# for generating simulation objects based on a config file

type GenConfig
	outputPath::String
	mode::String # "all" or "calls"
	
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
	targetResponseTime::Float
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
	transferDistrRng::DistrRng
	transferDurationDistrRng::DistrRng
	
	# misc RNGs
	ambStationRng::AbstractRNG
	callLocRng::AbstractRNG
	hospitalLocRng::AbstractRNG
	stationLocRng::AbstractRNG
	
	travelTimeFactorDistrRng::DistrRng
	
	GenConfig() = new("", "",
		"", "", "", "", "", "", "", "", "",
		nullIndex, nullIndex, nullIndex, nullIndex,
		nullIndex, nullIndex,
		Map(), 1e-6,
		nullTime, nullTime, nullTime, nullIndex, [],
		"", false, nullIndex, nullIndex)
end

function readGenConfig(genConfigFilename::String)
	# read gen config xml file
	rootElt = xmlFileRoot(genConfigFilename)
	@assert(name(rootElt) == "genConfig", string("xml root has incorrect name: ", name(rootElt)))
	
	genConfig = GenConfig()
	
	genConfig.outputPath = abspath(eltContentInterpVal(rootElt, "outputPath"))
	genConfig.mode = eltContent(rootElt, "mode")
	
	# output filenames
	simFilesElt = findElt(rootElt, "simFiles")
	simFilePath(filename::String) = joinpath(genConfig.outputPath, eltContent(simFilesElt, filename))
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
		seed = (seedAttr == nothing ? nullIndex : eval(parse(seedAttr)))
		return DistrRng(distr; seed = seed)
	end
	genConfig.interarrivalTimeDistrRng = callDistrsEltContent("interarrivalTime")
	genConfig.priorityDistrRng = callDistrsEltContent("priority")
	genConfig.dispatchDelayDistrRng = callDistrsEltContent("dispatchDelay")
	genConfig.onSceneDurationDistrRng = callDistrsEltContent("onSceneDuration")
	genConfig.transferDistrRng = callDistrsEltContent("transfer")
	genConfig.transferDurationDistrRng = callDistrsEltContent("transferDuration")
	
	# number of ambulances, calls, hospitals, stations
	genConfig.numAmbs = eltContentVal(simElt, "numAmbs")
	genConfig.numCalls = eltContentVal(simElt, "numCalls")
	genConfig.numHospitals = eltContentVal(simElt, "numHospitals")
	genConfig.numStations = eltContentVal(simElt, "numStations")
	
	# number of nodes in x and y direction for grid shaped graph
	graphElt = findElt(simElt, "graph")
	genConfig.xNodes = eltContentVal(graphElt, "xNodes")
	genConfig.yNodes = eltContentVal(graphElt, "yNodes")
	
	# misc values
	genConfig.startTime = eltContentVal(simElt, "startTime")
	@assert(genConfig.startTime >= 0)
	genConfig.targetResponseTime = eltContentVal(simElt, "targetResponseTime")
	genConfig.offRoadSpeed = eltContentVal(simElt, "offRoadSpeed") # km / day
	genConfig.stationCapacity = eltContentVal(simElt, "stationCapacity")
	travelModeSpeedsElt = findElt(simElt, "travelModeSpeeds")
	genConfig.travelModeSpeeds = (travelModeSpeedsElt == nothing ? [1.5 * genConfig.offRoadSpeed] : eltContentVal(travelModeSpeedsElt))
	
	# call gen parameters
	# call density raster
	callDensityRasterElt = findElt(simElt, "callDensityRaster")
	genConfig.callDensityRasterFilename = abspath(eltContentInterpVal(callDensityRasterElt, "filename"))
	genConfig.cropRaster = eltContentVal(callDensityRasterElt, "cropRaster")
	# seeds
	function callRasterSeedVal(seedName::String)
		seedAttr = attribute(callDensityRasterElt, seedName)
		return seedAttr == nothing ? nullIndex : eval(parse(seedAttr))
	end
	genConfig.callRasterCellSeed = callRasterSeedVal("cellSeed")
	genConfig.callRasterCellLocSeed = callRasterSeedVal("cellLocSeed")
	
	# some defaults - should move to config file sometime
	genConfig.ambStationRng = MersenneTwister(0)
	genConfig.callLocRng = MersenneTwister(1)
	genConfig.hospitalLocRng = MersenneTwister(4)
	genConfig.stationLocRng = MersenneTwister(5)
	genConfig.travelTimeFactorDistrRng = DistrRng(Uniform(1.0, 1.1); seed = 99)
	
	return genConfig
end

function runGenConfig(genConfigFilename::String; overwriteOutputPath::Bool = false, doPrint::Bool = true)
	genConfig = readGenConfig(genConfigFilename)
	
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
		writePrioritiesFile(genConfig.prioritiesFilename, repmat([genConfig.targetResponseTime],3))
		writeStationsFile(genConfig.stationsFilename, stations)
		writeTravelFile(genConfig.travelFilename, travel)
	elseif genConfig.mode == "calls"
		runGenConfigCalls(genConfig, doPrint = doPrint)
	else
		error("Unrecognised generation mode")
	end
end

# generate calls, and generate locations using call density raster file (file name stored in genConfig)
# raster may be cropped to be within map borders
function runGenConfigCalls(genConfig::GenConfig; doPrint::Bool = true)
	
	calls = makeCalls(genConfig) # will later overwrite location field values for each call
	numCalls = length(calls) # shorthand
	
	# read call density raster file
	doPrint && println("Reading raster file: ", genConfig.callDensityRasterFilename)
	raster = readRasterFile(genConfig.callDensityRasterFilename)
	
	if genConfig.cropRaster
		# crop raster to be within genConfig.map
		doPrint && println("Before cropping:")
		printRasterSize(raster)
		mapTrimmed = trimmedMap(genConfig.map, genConfig.mapTrim)
		if cropRaster!(raster, mapTrimmed)
			doPrint && println("After cropping:")
			printRasterSize(raster)
		end
	end
	
	rasterSampler = RasterSampler(raster, genConfig.callRasterCellSeed, genConfig.callRasterCellLocSeed)
	randLocations = rasterRandLocations(rasterSampler, numCalls)
	for i = 1:numCalls
		calls[i].location = randLocations[i]
	end
	
	doPrint && println("Saving calls file to: ", genConfig.callsFilename)
	writeCallsFile(genConfig.callsFilename, genConfig.startTime, calls)
end

function makeAmbs(genConfig::GenConfig)
	ambulances = Vector{Ambulance}(genConfig.numAmbs)
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
	
	arcs = Vector{Arc}(graph.ne)
	travelTimes = Array{Float,2}(numTravelModes,length(arcs))
	
	i = 1
	for edge in LightGraphs.edges(graph)
		arcs[i] = Arc()
		arcs[i].index = i
		arcs[i].fromNodeIndex = edge.src
		arcs[i].toNodeIndex = edge.dst
		i = i + 1
	end
	
	for i = 1:numTravelModes, arc in arcs
		dist = normDist(genConfig.map, nodes[arc.fromNodeIndex].location, nodes[arc.toNodeIndex].location)
		travelTimes[i,arc.index] = dist / genConfig.travelModeSpeeds[i] * rand(genConfig.travelTimeFactorDistrRng)
	end
	
	return arcs, travelTimes
end

# make calls that are spatially randomly uniform in map
function makeCalls(genConfig::GenConfig)
	calls = Vector{Call}(genConfig.numCalls)
	
	currentTime = genConfig.startTime
	# first call will arrive at genConfig.startTime + rand(genConfig.interarrivalTimeDistrRng)
	for i = 1:genConfig.numCalls
		currentTime += rand(genConfig.interarrivalTimeDistrRng) # apply time step
		
		calls[i] = Call()
		calls[i].index = i
		calls[i].priority = Priority(rand(genConfig.priorityDistrRng))
		calls[i].location = randLocation(genConfig.map; trim = genConfig.mapTrim, rng = genConfig.callLocRng)
		calls[i].arrivalTime = currentTime
		calls[i].dispatchDelay = rand(genConfig.dispatchDelayDistrRng)
		calls[i].onSceneDuration = rand(genConfig.onSceneDurationDistrRng)
		calls[i].transfer = (rand(genConfig.transferDistrRng) == 1)
		calls[i].hospitalIndex = nullIndex
		calls[i].transferDuration = rand(genConfig.transferDurationDistrRng)
	end
	
	return calls
end

function makeHospitals(genConfig::GenConfig)
	hospitals = Vector{Hospital}(genConfig.numHospitals)
	for i = 1:genConfig.numHospitals
		hospitals[i] = Hospital()
		hospitals[i].index = i
		hospitals[i].location = randLocation(genConfig.map; trim = genConfig.mapTrim, rng = genConfig.hospitalLocRng)
	end
	return hospitals
end

function makeNodes(genConfig::GenConfig, graph::LightGraphs.Graph)
	nodes = Vector{Node}(LightGraphs.nv(graph)) # should have length xNodes*yNodes
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
	stations = Vector{Station}(genConfig.numStations)
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
	# all three travel priorities should be used
	travel.numSets = 1
	travel.modeLookup = repmat([1],1,3) # modeLookup[setIndex, priority] gives travelModeIndex
	
	travel.setsStartTimes = [genConfig.startTime]
	travel.setsTimeOrder = [1]
	
	return travel
end
