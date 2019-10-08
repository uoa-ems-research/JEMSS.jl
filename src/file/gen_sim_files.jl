##########################################################################
# Copyright 2017 Samuel Ridler.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
##########################################################################

# for generating simulation objects based on a config file

mutable struct GenConfig
	inputPath::String
	outputPath::String
	mode::String # "all" or "calls"
	numCallsFiles::Int
	
	# output file names
	ambsFilename::String
	arcsFilename::String
	callsFilename::String
	hospitalsFilename::String
	mapFilename::String
	nodesFilename::String
	prioritiesFilename::String
	stationsFilename::String
	travelFilename::String
	
	# counts
	numAmbs::Int
	numCalls::Int
	numHospitals::Int
	numStations::Int
	
	# graph
	xNodes::Int # number of nodes in x direction
	yNodes::Int # number of nodes in y direction
	
	# map
	map::Map
	mapTrim::Float # fraction of map border to trim, to make sure objects generated on map are inside borders
	
	# misc
	startTime::Float
	maxCallArrivalTime::Float
	minLastCallArrivalTime::Float
	targetResponseDuration::Float
	offRoadSpeed::Float
	stationCapacity::Int
	travelModeSpeeds::Vector{Float}
	
	# call density raster
	callDensityRasterFilename::String
	cropRaster::Bool
	callRasterCellSeed::Int # seed for rng, will generate raster cell index
	callRasterCellLocSeed::Int # seed for rng, will generate location within raster cell
	
	# call related distributions and random number generators
	interarrivalTimeDistrRng::DistrRng
	priorityDistrRng::DistrRng
	dispatchDelayDistrRng::DistrRng
	onSceneDurationDistrRng::DistrRng
	transportDistrRng::DistrRng
	handoverDurationDistrRng::DistrRng
	
	# misc RNGs
	ambStationRng::AbstractRNG
	callLocRng::AbstractRNG
	hospitalLocRng::AbstractRNG
	stationLocRng::AbstractRNG
	
	# undefined
	travelTimeFactorDistrRng::DistrRng
	callRasterSampler::RasterSampler
	
	GenConfig() = new("", "", "", 1,
		"", "", "", "", "", "", "", "", "",
		nullIndex, nullIndex, nullIndex, nullIndex,
		nullIndex, nullIndex,
		Map(), 1e-6,
		nullTime, nullTime, nullTime, nullTime, nullTime, nullIndex, [],
		"", false, nullIndex, nullIndex)
end

function readGenConfig(genConfigFilename::String)
	# read gen config xml file
	genConfigFilename = genConfigFilename |> interpolateString |> abspath
	global genConfigFileDir = dirname(genConfigFilename)
	rootElt = xmlFileRoot(genConfigFilename)
	@assert(xName(rootElt) == "genConfig", string("xml root has incorrect name: ", xName(rootElt)))
	
	genConfig = GenConfig()
	
	genConfig.inputPath = containsElt(rootElt, "inputPath") ? joinPathIfNotAbs(genConfigFileDir, eltContentInterpVal(rootElt, "inputPath")) : genConfigFileDir
	genConfig.outputPath = joinPathIfNotAbs(genConfigFileDir, eltContentInterpVal(rootElt, "outputPath")) # output path can be absolute, or relative to genConfigFileDir
	genConfig.mode = eltContent(rootElt, "mode")
	genConfig.numCallsFiles = containsElt(rootElt, "numCallsFiles") ? eltContentVal(rootElt, "numCallsFiles") : 1
	@assert(genConfig.numCallsFiles >= 1)
	
	# output filenames
	simFilesElt = findElt(rootElt, "simFiles")
	simFilePath(filename::String) = joinPathIfNotAbs(genConfig.outputPath, eltContentInterpVal(simFilesElt, filename)) # filename can be absolute, or relative to output path
	genConfig.ambsFilename = simFilePath("ambulances")
	genConfig.arcsFilename = simFilePath("arcs")
	genConfig.callsFilename = simFilePath("calls")
	genConfig.hospitalsFilename = simFilePath("hospitals")
	genConfig.mapFilename = simFilePath("map")
	genConfig.nodesFilename = simFilePath("nodes")
	genConfig.prioritiesFilename = simFilePath("priorities")
	genConfig.stationsFilename = simFilePath("stations")
	genConfig.travelFilename = simFilePath("travel")
	
	# read sim parameters
	simElt = findElt(rootElt, "sim")
	
	# create map
	# map is needed before generating random locations
	mapElt = findElt(simElt, "map")
	map = Map()
	map.xMin = eltContentVal(mapElt, "xMin")
	map.xMax = eltContentVal(mapElt, "xMax")
	map.xScale = eltContentVal(mapElt, "xScale")
	map.xRange = map.xMax - map.xMin
	map.yMin = eltContentVal(mapElt, "yMin")
	map.yMax = eltContentVal(mapElt, "yMax")
	map.yScale = eltContentVal(mapElt, "yScale")
	map.yRange = map.yMax - map.yMin
	@assert(map.xRange > 0 && map.yRange > 0)
	genConfig.map = map
	
	# call distributions and random number generators
	callDistrsElt = findElt(simElt, "callDistributions")
	function callDistrsEltContent(distrName::String)
		distrElt = findElt(callDistrsElt, distrName)
		distr = eltContentVal(distrElt)
		seedAttr = attribute(distrElt, "seed")
		seed = (seedAttr == nothing ? nullIndex : eval(Meta.parse(seedAttr)))
		return DistrRng(distr; seed = seed)
	end
	genConfig.interarrivalTimeDistrRng = callDistrsEltContent("interarrivalTime")
	genConfig.priorityDistrRng = callDistrsEltContent("priority")
	genConfig.dispatchDelayDistrRng = callDistrsEltContent("dispatchDelay")
	genConfig.onSceneDurationDistrRng = callDistrsEltContent("onSceneDuration")
	if containsElt(callDistrsElt, "transport") && containsElt(callDistrsElt, "handoverDuration")
		genConfig.transportDistrRng = callDistrsEltContent("transport")
		genConfig.handoverDurationDistrRng = callDistrsEltContent("handoverDuration")
	else # compat
		@assert(containsElt(callDistrsElt, "transfer") && containsElt(callDistrsElt, "transferDuration"))
		@warn("`transfer` and `transferDuration` elements for call distributions in gen config are deprecated, use `transport` and `handoverDuration` instead.")
		genConfig.transportDistrRng = callDistrsEltContent("transfer")
		genConfig.handoverDurationDistrRng = callDistrsEltContent("transferDuration")
	end
	
	# number of ambulances, calls, hospitals, stations
	genConfig.numAmbs = eltContentVal(simElt, "numAmbs")
	genConfig.numCalls = containsElt(simElt, "numCalls") ? eltContentVal(simElt, "numCalls") : nullIndex # can alternatively specify maxCallArrivalTime or minLastCallArrivalTime
	genConfig.numHospitals = eltContentVal(simElt, "numHospitals")
	genConfig.numStations = eltContentVal(simElt, "numStations")
	
	# number of nodes in x and y direction for grid shaped graph
	graphElt = findElt(simElt, "graph")
	genConfig.xNodes = eltContentVal(graphElt, "xNodes")
	genConfig.yNodes = eltContentVal(graphElt, "yNodes")
	
	# misc values
	genConfig.startTime = eltContentVal(simElt, "startTime")
	@assert(genConfig.startTime >= 0)
	genConfig.maxCallArrivalTime = containsElt(simElt, "maxCallArrivalTime") ? eltContentVal(simElt, "maxCallArrivalTime") : nullTime
	genConfig.minLastCallArrivalTime = containsElt(simElt, "minLastCallArrivalTime") ? eltContentVal(simElt, "minLastCallArrivalTime") : nullTime
	@assert((genConfig.numCalls != nullIndex) + (genConfig.maxCallArrivalTime != nullTime) + (genConfig.minLastCallArrivalTime != nullTime) == 1, "Need exactly one of these values: numCalls, maxCallArrivalTime, minLastCallArrivalTime.")
	@assert(genConfig.startTime <= genConfig.maxCallArrivalTime || genConfig.maxCallArrivalTime == nullTime)
	@assert(genConfig.startTime <= genConfig.minLastCallArrivalTime || genConfig.minLastCallArrivalTime == nullTime)
	targetResponseDurationString = containsElt(simElt, "targetResponseDuration") ? "targetResponseDuration" : "targetResponseTime" # compat for "targetResponseTime"
	genConfig.targetResponseDuration = eltContentVal(simElt, targetResponseDurationString)
	genConfig.offRoadSpeed = eltContentVal(simElt, "offRoadSpeed") # km / day
	genConfig.stationCapacity = eltContentVal(simElt, "stationCapacity")
	travelModeSpeedsElt = findElt(simElt, "travelModeSpeeds")
	genConfig.travelModeSpeeds = (travelModeSpeedsElt == nothing ? [1.5 * genConfig.offRoadSpeed] : eltContentVal(travelModeSpeedsElt))
	
	# call gen parameters
	# call density raster
	callDensityRasterElt = findElt(simElt, "callDensityRaster")
	genConfig.callDensityRasterFilename = joinPathIfNotAbs(genConfig.inputPath, eltContentInterpVal(callDensityRasterElt, "filename"))
	genConfig.cropRaster = eltContentVal(callDensityRasterElt, "cropRaster")
	# seeds
	function callRasterSeedVal(seedName::String)
		seedAttr = attribute(callDensityRasterElt, seedName)
		return seedAttr == nothing ? nullIndex : eval(Meta.parse(seedAttr))
	end
	genConfig.callRasterCellSeed = callRasterSeedVal("cellSeed")
	genConfig.callRasterCellLocSeed = callRasterSeedVal("cellLocSeed")
	if isfile(genConfig.callDensityRasterFilename)
		raster = readRasterFile(genConfig.callDensityRasterFilename)
		if genConfig.cropRaster
			# crop raster to be within genConfig.map
			mapTrimmed = trimmedMap(genConfig.map, genConfig.mapTrim)
			cropRaster!(raster, mapTrimmed)
		end
		genConfig.callRasterSampler = RasterSampler(raster, genConfig.callRasterCellSeed, genConfig.callRasterCellLocSeed)
	end
	
	# some defaults - should move to config file sometime
	genConfig.ambStationRng = MersenneTwister(0)
	genConfig.callLocRng = MersenneTwister(1)
	genConfig.hospitalLocRng = MersenneTwister(4)
	genConfig.stationLocRng = MersenneTwister(5)
	genConfig.travelTimeFactorDistrRng = DistrRng(Uniform(1.0, 1.1); seed = 99)
	
	return genConfig
end

function runGenConfig(genConfig::GenConfig; overwriteOutputPath::Bool = false, doPrint::Bool = true)
	if isdir(genConfig.outputPath) && !overwriteOutputPath
		doPrint && println("Output path already exists: ", genConfig.outputPath)
		doPrint && print("Delete folder contents and continue anyway? (y = yes): ")
		response = chomp(readline())
		if response != "y"
			doPrint && println("stopping")
			return
		else
			overwriteOutputPath = true
		end
	end
	if isdir(genConfig.outputPath) && overwriteOutputPath
		doPrint && println("Deleting folder contents: ", genConfig.outputPath)
		rm(genConfig.outputPath; recursive=true)
	end
	if !isdir(genConfig.outputPath)
		mkdir(genConfig.outputPath)
	end
	
	doPrint && println("Generation mode: ", genConfig.mode)
	if genConfig.mode == "all"
		# make all
		ambulances = makeAmbs(genConfig)
		calls = makeCalls(genConfig)
		hospitals = makeHospitals(genConfig)
		stations = makeStations(genConfig)
		travel = makeTravel(genConfig)
		
		graph = LightGraphs.Grid([genConfig.xNodes, genConfig.yNodes]) # grid shaped graph
		nodes = makeNodes(genConfig, graph)
		(arcs, travelTimes) = makeArcs(genConfig, graph, nodes)
		
		# save all
		doPrint && println("Saving output to: ", genConfig.outputPath)
		writeAmbsFile(genConfig.ambsFilename, ambulances)
		writeArcsFile(genConfig.arcsFilename, arcs, travelTimes, "undirected")
		writeCallsFile(genConfig.callsFilename, genConfig.startTime, calls)
		writeHospitalsFile(genConfig.hospitalsFilename, hospitals)
		writeMapFile(genConfig.mapFilename, genConfig.map)
		writeNodesFile(genConfig.nodesFilename, nodes)
		writePrioritiesFile(genConfig.prioritiesFilename, repeat([genConfig.targetResponseDuration],numPriorities))
		writeStationsFile(genConfig.stationsFilename, stations)
		writeTravelFile(genConfig.travelFilename, travel)
	elseif genConfig.mode == "calls"
		runGenConfigCalls(genConfig, doPrint = doPrint)
	else
		error("Unrecognised generation mode")
	end
	
	return genConfig
end
function runGenConfig(genConfigFilename::String; overwriteOutputPath::Bool = false, doPrint::Bool = true)
	genConfig = readGenConfig(genConfigFilename)
	return runGenConfig(genConfig, overwriteOutputPath = overwriteOutputPath, doPrint = doPrint)
end

# generate calls, and generate locations using call density raster file (file name stored in genConfig)
function runGenConfigCalls(genConfig::GenConfig; doPrint::Bool = true, writeFile::Bool = true)
	numFiles = genConfig.numCallsFiles # shorthand
	@assert(numFiles >= 1)
	doPrint && println("Generating $numFiles call set(s).")
	callSets = [makeCalls(genConfig) for i = 1:numFiles]
	
	if writeFile
		if numFiles == 1
			calls = callSets[1]
			filename = genConfig.callsFilename
			doPrint && println("Saving calls file to: ", filename)
			writeCallsFile(filename, genConfig.startTime, calls)
		else
			@assert(numFiles > 1)
			doPrint && println("Saving calls files to: ", genConfig.outputPath)
			for j = 1:numFiles
				filename = joinpath(genConfig.outputPath, "calls_$j.csv")
				writeCallsFile(filename, genConfig.startTime, callSets[j])
			end
		end
	end
	
	return callSets
end

function makeAmbs(genConfig::GenConfig)
	ambulances = Vector{Ambulance}(undef, genConfig.numAmbs)
	for i = 1:genConfig.numAmbs
		ambulances[i] = Ambulance()
		ambulances[i].index = i
		ambulances[i].stationIndex = rand(genConfig.ambStationRng, 1:genConfig.numStations)
		ambulances[i].class = als
	end
	return ambulances
end

function makeArcs(genConfig::GenConfig, graph::LightGraphs.Graph, nodes::Vector{Node})
	numTravelModes = length(genConfig.travelModeSpeeds)
	@assert(numTravelModes >= 1)
	
	arcs = Vector{Arc}(undef, graph.ne)
	travelTimes = Array{Float,2}(undef,numTravelModes,length(arcs))
	
	i = 1
	for edge in LightGraphs.edges(graph)
		arcs[i] = Arc()
		arcs[i].index = i
		arcs[i].fromNodeIndex = edge.src
		arcs[i].toNodeIndex = edge.dst
		i = i + 1
	end
	
	for i = 1:numTravelModes, arc in arcs
		arc.distance = normDist(genConfig.map, nodes[arc.fromNodeIndex].location, nodes[arc.toNodeIndex].location)
		travelTimes[i,arc.index] = arc.distance / genConfig.travelModeSpeeds[i] * rand(genConfig.travelTimeFactorDistrRng)
	end
	
	return arcs, travelTimes
end

# make calls that are spatially randomly uniform in map, or distributed according to raster
function makeCalls(genConfig::GenConfig; rasterSampler::Union{RasterSampler,Nothing} = nothing)
	if (rasterSampler == nothing && isdefined(genConfig, :callRasterSampler)) rasterSampler = genConfig.callRasterSampler end
	calls = Call[]
	currentTime = genConfig.startTime
	while length(calls) < genConfig.numCalls || genConfig.numCalls == nullIndex
		currentTime += rand(genConfig.interarrivalTimeDistrRng) # apply time step
		if (currentTime > genConfig.maxCallArrivalTime && genConfig.maxCallArrivalTime != nullTime) break end
		
		call = Call()
		call.index = length(calls) + 1
		call.priority = Priority(rand(genConfig.priorityDistrRng))
		call.arrivalTime = currentTime
		call.dispatchDelay = rand(genConfig.dispatchDelayDistrRng)
		call.onSceneDuration = rand(genConfig.onSceneDurationDistrRng)
		call.transport = (rand(genConfig.transportDistrRng) == 1)
		call.hospitalIndex = nullIndex
		call.handoverDuration = rand(genConfig.handoverDurationDistrRng)
		if !call.transport call.handoverDuration = 0.0 end
		if rasterSampler == nothing
			call.location = randLocation(genConfig.map; trim = genConfig.mapTrim, rng = genConfig.callLocRng)
		else
			call.location = rasterRandLocations(rasterSampler, 1)[1]
		end
		
		push!(calls, call)
		
		if (currentTime >= genConfig.minLastCallArrivalTime && genConfig.minLastCallArrivalTime != nullTime) break end
	end
	
	return calls
end

function makeHospitals(genConfig::GenConfig)
	hospitals = Vector{Hospital}(undef, genConfig.numHospitals)
	for i = 1:genConfig.numHospitals
		hospitals[i] = Hospital()
		hospitals[i].index = i
		hospitals[i].location = randLocation(genConfig.map; trim = genConfig.mapTrim, rng = genConfig.hospitalLocRng)
	end
	return hospitals
end

function makeNodes(genConfig::GenConfig, graph::LightGraphs.Graph)
	nodes = Vector{Node}(undef, LightGraphs.nv(graph)) # should have length xNodes*yNodes
	map = genConfig.map # shorthand
	k = 1
	for j = 1:genConfig.yNodes, i = 1:genConfig.xNodes
		nodes[k] = Node()
		nodes[k].index = k
		nodes[k].location.x = map.xRange * ((i-0.5)/genConfig.xNodes) + map.xMin
		nodes[k].location.y = map.yRange * ((j-0.5)/genConfig.yNodes) + map.yMin
		k += 1
	end
	
	return nodes
end

function makeStations(genConfig::GenConfig)
	stations = Vector{Station}(undef, genConfig.numStations)
	for i = 1:genConfig.numStations
		stations[i] = Station()
		stations[i].index = i
		stations[i].location = randLocation(genConfig.map; trim = genConfig.mapTrim, rng = genConfig.stationLocRng)
		stations[i].capacity = genConfig.stationCapacity
	end
	return stations
end

function makeTravel(genConfig::GenConfig)
	travel = Travel() # note: will only partially fill travel fields, just enough for writeTravelFile()
	
	# create single travel mode, add to travel
	travelMode = TravelMode()
	travelMode.index = 1
	travelMode.offRoadSpeed = genConfig.offRoadSpeed
	travel.modes = [travelMode]
	travel.numModes = length(travel.modes)
	
	# use single travel mode in single travel mode set
	travel.numSets = 1
	travel.modeLookup = repeat([1],1,numPriorities) # modeLookup[setIndex, priority] gives travelModeIndex
	
	travel.setsStartTimes = [genConfig.startTime]
	travel.setsTimeOrder = [1]
	
	return travel
end
