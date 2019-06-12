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

# run configuration xml file
function runConfig(configFilename::String)
	# read input files, initialise simulation
	sim = initSim(configFilename; allowResim = true, createBackup = false, allowWriteOutput = true)
	
	# open output files
	openOutputFiles!(sim)
	
	simulate!(sim)
	
	# print statistics
	printSimStats(sim)
	
	# save statistics
	writeStatsFiles!(sim)
	
	# close output files
	closeOutputFiles!(sim)
end

# initialise simulation from input files
function initSim(configFilename::String;
	allowResim::Bool = false, createBackup::Bool = true, allowWriteOutput::Bool = false, doPrint::Bool = true)
	
	# read sim config xml file
	configFilename = configFilename |> interpolateString |> abspath
	configFileDir = splitdir(configFilename)[1]
	rootElt = xmlFileRoot(configFilename)
	@assert(xName(rootElt) == "simConfig", string("xml root has incorrect name: ", xName(rootElt)))
	
	# for progress messages:
	t = [0.0]
	initMessage(t, msg) = doPrint && (t[1] = time(); print(msg))
	initTime(t) = doPrint && println(": ", round(time() - t[1], digits = 2), " seconds")
	
	##################
	# sim config
	
	initMessage(t, "reading config file data")
	
	sim = Simulation()
	sim.configRootElt = rootElt
	
	# input
	sim.inputPath = joinPathIfNotAbs(configFileDir, eltContentInterpVal(rootElt, "inputPath")) # input path can be absolute, or relative to configFilename
	simFilesElt = findElt(rootElt, "simFiles")
	inputFiles = childrenNodeNames(simFilesElt)
	sim.inputFiles = Dict{String,File}()
	for inputFile in inputFiles
		file = File()
		file.path = joinPathIfNotAbs(sim.inputPath, eltContentInterpVal(simFilesElt, inputFile))
		file.name = splitdir(file.path)[2]
		if inputFile != "rNetTravels" # do not need checksum of rNetTravels file
			file.checksum = fileChecksum(file.path)
		end
		sim.inputFiles[inputFile] = file
	end
	
	# output
	sim.writeOutput = allowWriteOutput && eltContentVal(rootElt, "writeOutput")
	sim.outputPath = joinPathIfNotAbs(configFileDir, eltContentInterpVal(rootElt, "outputPath")) # output path can be absolute, or relative to configFilename
	outputFilesElt = findElt(rootElt, "outputFiles")
	outputFiles = childrenNodeNames(outputFilesElt)
	sim.outputFiles = Dict{String,File}()
	for outputFile in outputFiles
		file = File()
		file.path = joinPathIfNotAbs(sim.outputPath, eltContentInterpVal(outputFilesElt, outputFile))
		file.name = splitdir(file.path)[2]
		sim.outputFiles[outputFile] = file
	end
	
	initTime(t)
	
	##################
	# read simulation input files
	
	initMessage(t, "reading input data")
	
	simFilePath(name::String) = sim.inputFiles[name].path
	
	# read sim data
	sim.ambulances = readAmbsFile(simFilePath("ambulances"))
	(sim.calls, sim.startTime) = readCallsFile(simFilePath("calls"))
	sim.time = sim.startTime
	sim.hospitals = readHospitalsFile(simFilePath("hospitals"))
	sim.stations = readStationsFile(simFilePath("stations"))
	
	sim.numAmbs = length(sim.ambulances)
	sim.numCalls = length(sim.calls)
	sim.numHospitals = length(sim.hospitals)
	sim.numStations = length(sim.stations)
	
	# read network data
	sim.net = Network()
	net = sim.net # shorthand
	fGraph = net.fGraph # shorthand
	fGraph.nodes = readNodesFile(simFilePath("nodes"))
	(fGraph.arcs, arcTravelTimes) = readArcsFile(simFilePath("arcs"))
	
	# read rNetTravels from file, if saved
	rNetTravelsLoaded = NetTravel[]
	rNetTravelsFilename = ""
	if haskey(sim.inputFiles, "rNetTravels")
		rNetTravelsFilename = simFilePath("rNetTravels")
		if isfile(rNetTravelsFilename)
			rNetTravelsLoaded = readRNetTravelsFile(rNetTravelsFilename)
		elseif !isdir(dirname(rNetTravelsFilename)) || splitdir(rNetTravelsFilename)[2] == ""
			# rNetTravelsFilename is invalid
			rNetTravelsFilename = ""
		else
			# save net.rNetTravels to file once calculated
		end
	end
	
	# read misc
	sim.map = readMapFile(simFilePath("map"))
	map = sim.map # shorthand
	(sim.targetResponseDurations, sim.responseTravelPriorities) = readPrioritiesFile(simFilePath("priorities"))
	sim.travel = readTravelFile(simFilePath("travel"))
	# "demand" file can be slow to read, will not read here but read elsewhere when needed
	if haskey(sim.inputFiles, "demandCoverage")
		sim.demandCoverage = readDemandCoverageFile(simFilePath("demandCoverage"))
	end
	
	initTime(t)
	
	##################
	# network
	
	initMessage(t, "initialising fGraph")
	initGraph!(fGraph)
	initTime(t)
	
	initMessage(t, "checking fGraph")
	checkGraph(fGraph, map)
	initTime(t)
	
	initMessage(t, "initialising fNetTravels")
	initFNetTravels!(net, arcTravelTimes)
	initTime(t)
	
	initMessage(t, "creating rGraph from fGraph")
	createRGraphFromFGraph!(net)
	initTime(t)
	doPrint && println("fNodes: ", length(net.fGraph.nodes), ", rNodes: ", length(net.rGraph.nodes))
	
	initMessage(t, "checking rGraph")
	checkGraph(net.rGraph, map)
	initTime(t)
	
	if rNetTravelsLoaded != []
		doPrint && println("using data from rNetTravels file")
		try
			initMessage(t, "creating rNetTravels from fNetTravels")
			createRNetTravelsFromFNetTravels!(net; rNetTravelsLoaded = rNetTravelsLoaded)
			initTime(t)
		catch
			doPrint && println()
			@warn("failed to use data from rNetTravels file")
			rNetTravelsLoaded = []
			rNetTravelsFilename = ""
		end
	end
	if rNetTravelsLoaded == []
		initMessage(t, "creating rNetTravels from fNetTravels, and shortest paths")
		createRNetTravelsFromFNetTravels!(net)
		initTime(t)
		if rNetTravelsFilename != ""
			initMessage(t, "saving rNetTravels to file")
			writeRNetTravelsFile(rNetTravelsFilename, net.rNetTravels)
			initTime(t)
		end
	end
	
	##################
	# travel
	
	initMessage(t, "initialising travel")
	
	travel = sim.travel # shorthand
	@assert(travel.setsStartTimes[1] <= sim.startTime)
	@assert(length(net.fNetTravels) == travel.numModes)
	for travelMode in travel.modes
		travelMode.fNetTravel = net.fNetTravels[travelMode.index]
		travelMode.rNetTravel = net.rNetTravels[travelMode.index]
	end
	
	initTime(t)
	
	##################
	# grid
	
	initMessage(t, "placing nodes in grid")
	
	# hard-coded grid size
	# grid rects will be roughly square, with one node per square on average
	n = length(fGraph.nodes)
	xDist = map.xRange * map.xScale
	yDist = map.yRange * map.yScale
	nx = Int(ceil(sqrt(n * xDist / yDist)))
	ny = Int(ceil(sqrt(n * yDist / xDist)))
	
	sim.grid = Grid(map, nx, ny)
	grid = sim.grid # shorthand
	gridPlaceNodes!(map, grid, fGraph.nodes)
	initTime(t)
	
	doPrint && println("nodes: ", length(fGraph.nodes), ", grid size: ", nx, " x ", ny)
	
	##################
	# sim - ambulances, calls, hospitals, stations...
	
	initMessage(t, "adding ambulances, calls, etc")
	
	# for each call, hospital, and station, find nearest node
	for c in sim.calls
		(c.nearestNodeIndex, c.nearestNodeDist) = findNearestNodeInGrid(map, grid, fGraph.nodes, c.location)
	end
	for h in sim.hospitals
		(h.nearestNodeIndex, h.nearestNodeDist) = findNearestNodeInGrid(map, grid, fGraph.nodes, h.location)
	end
	for s in sim.stations
		(s.nearestNodeIndex, s.nearestNodeDist) = findNearestNodeInGrid(map, grid, fGraph.nodes, s.location)
	end
	
	# create event list
	# try to add events to eventList in reverse time order, to reduce sorting required
	sim.eventList = Vector{Event}()
	
	# add first call to event list
	addEvent!(sim.eventList, sim.calls[1])
	
	# create ambulance wake up events
	for a in sim.ambulances
		initAmbulance!(sim, a)
		# currently, this sets ambulances to wake up at start of sim, since wake up and sleep events are not in ambulances file yet
	end
	
	initTime(t)
	
	initMessage(t, "storing times between fNodes and common locations")
	
	# for each station, find time to each node in fGraph for each travel mode, and vice versa (node to station)
	# for each node in fGraph and each travel mode, find nearest hospital
	# requires deterministic and static travel times
	
	commonFNodes = sort(unique(vcat([h.nearestNodeIndex for h in sim.hospitals], [s.nearestNodeIndex for s in sim.stations])))
	setCommonFNodes!(net, commonFNodes)
	
	# find the nearest hospital to travel to from each node in fGraph
	numFNodes = length(fGraph.nodes) # shorthand
	for fNetTravel in net.fNetTravels
		fNetTravel.fNodeNearestHospitalIndex = Vector{Int}(undef, numFNodes)
		travelModeIndex = fNetTravel.modeIndex # shorthand
		travelMode = travel.modes[travelModeIndex] # shorthand
		for node in fGraph.nodes
			# find nearest hospital to node
			minTime = Inf
			nearestHospitalIndex = nullIndex
			for hospital in sim.hospitals
				travelTime = shortestPathTravelTime(net, travelModeIndex, node.index, hospital.nearestNodeIndex)
				travelTime += offRoadTravelTime(travelMode, hospital.nearestNodeDist)
				if travelTime < minTime
					minTime = travelTime
					nearestHospitalIndex = hospital.index
				end
			end
			fNetTravel.fNodeNearestHospitalIndex[node.index] = nearestHospitalIndex
		end
	end
	
	initTime(t)
	
	##################
	# decision logic
	
	decisionElt = findElt(rootElt, "decision")
	sim.addCallToQueue! = eltContentVal(decisionElt, "callQueueing")
	sim.findAmbToDispatch! = eltContentVal(decisionElt, "dispatch")
	haskey(sim.inputFiles, "redispatch") && (sim.redispatch = readRedispatchFile(simFilePath("redispatch")))
	
	# move up
	mud = sim.moveUpData # shorthand
	moveUpElt = findElt(decisionElt, "moveUp")
	moveUpModuleName = eltContent(moveUpElt, "module")
	
	if moveUpModuleName == "none"
		mud.useMoveUp = false
		mud.moveUpModule = nullMoveUpModule
		doPrint && println("not using move up")
	else
		initMessage(t, "initialising move up")
		
		mud.useMoveUp = true
		
		if moveUpModuleName == "comp_table"
			mud.moveUpModule = compTableModule
			compTableElt = findElt(moveUpElt, "compTable")
			compTableFilename = joinPathIfNotAbs(sim.inputPath, eltContentInterpVal(compTableElt, "filename"))
			initCompTable!(sim, compTableFilename)
			
		elseif moveUpModuleName == "ddsm"
			mud.moveUpModule = ddsmModule
			ddsmElt = findElt(moveUpElt, "ddsm")
			initDdsm!(sim;
				coverFractionTargetT1 = eltContentVal(ddsmElt, "coverFractionTargetT1"),
				travelTimeCost = eltContentVal(ddsmElt, "travelTimeCost"),
				slackWeight = eltContentVal(ddsmElt, "slackWeight"),
				coverTimeDemandPriorities = eltContentVal(ddsmElt, "coverTimeDemandPriorities"),
				options = convert(Dict{Symbol,Any}, Dict(eltContentVal(ddsmElt, "options"))))
			
		elseif moveUpModuleName == "dmexclp"
			mud.moveUpModule = dmexclpModule
			dmexclpElt = findElt(moveUpElt, "dmexclp")
			initDmexclp!(sim; busyFraction = eltContentVal(dmexclpElt, "busyFraction"))
			
		elseif moveUpModuleName == "priority_list"
			mud.moveUpModule = priorityListModule
			priorityListElt = findElt(moveUpElt, "priorityList")
			priorityListFilename = joinPathIfNotAbs(sim.inputPath, eltContentInterpVal(priorityListElt, "filename"))
			initPriorityList!(sim, priorityListFilename)
			
		elseif moveUpModuleName == "zhang_ip"
			mud.moveUpModule = zhangIpModule
			zhangIpElt = findElt(moveUpElt, "zhangIp")
			zhangIpParamsFilename = joinPathIfNotAbs(sim.inputPath, eltContentInterpVal(zhangIpElt, "paramsFilename"))
			initZhangIp!(sim;
				paramsFilename = zhangIpParamsFilename)
			
		elseif moveUpModuleName == "temp0"
			mud.moveUpModule = temp0Module
			temp0Elt = findElt(moveUpElt, "temp0")
			initTemp0!(sim;
				busyFraction = eltContentVal(temp0Elt, "busyFraction"),
				travelTimeCost = eltContentVal(temp0Elt, "travelTimeCost"),
				maxIdleAmbTravelTime = eltContentVal(temp0Elt, "maxIdleAmbTravelTime"),
				maxNumNearestStations = eltContentVal(temp0Elt, "maxNumNearestStations"))
			
		elseif moveUpModuleName == "temp1"
			mud.moveUpModule = temp1Module
			temp1Elt = findElt(moveUpElt, "temp1")
			initTemp1!(sim;
				busyFraction = eltContentVal(temp1Elt, "busyFraction"),
				travelTimeCost = eltContentVal(temp1Elt, "travelTimeCost"),
				maxIdleAmbTravelTime = eltContentVal(temp1Elt, "maxIdleAmbTravelTime"),
				maxNumNearestStations = eltContentVal(temp1Elt, "maxNumNearestStations"))
			
		elseif moveUpModuleName == "temp2"
			mud.moveUpModule = temp2Module
			temp2Elt = findElt(moveUpElt, "temp2")
			initTemp2!(sim;
				busyFraction = eltContentVal(temp2Elt, "busyFraction"),
				travelTimeCost = eltContentVal(temp2Elt, "travelTimeCost"),
				maxIdleAmbTravelTime = eltContentVal(temp2Elt, "maxIdleAmbTravelTime"),
				maxNumNearestStations = eltContentVal(temp2Elt, "maxNumNearestStations"))
		else
			error("invalid move up module name given: ", moveUpModuleName)
		end
		
		initTime(t)
		
		doPrint && println("using move up module: ", moveUpModuleName)
	end
	
	##################
	# resimulation
	
	if allowResim
		sim.resim.use = eltContentVal(rootElt, "resim")
		if sim.resim.use
			initMessage(t, "")
			initResim!(sim)
			doPrint && print("initialised resimulation")
			initTime(t)
		end
	end
	
	##################
	# statistics
	
	initMessage(t, "initialising statistics")
	
	# temp
	@info("Todo: make sim stats capturing configurable.")
	stats = sim.stats # shorthand
	stats.doCapture = true
	stats.captureTimes = sim.startTime + 1 : 1 : sim.calls[end].arrivalTime + 1 # daily
	stats.nextCaptureTime = stats.captureTimes[1]
	
	initTime(t)
	doPrint && stats.doCapture && println("will capture statistics")
	
	##################
	# misc
	
	sim.initialised = true # at this point, the simulation could be run
	
	##################
	# sim backup
	
	if createBackup
		initMessage(t, "creating sim backup")
		backup!(sim) # for restarting sim
		initTime(t)
	end
	
	return sim
end

# initialise given ambulance
# sets ambulance as sleeping, creates wake up event
function initAmbulance!(sim::Simulation, ambulance::Ambulance;
	wakeUpTime::Float = nullTime)
	wakeUpTime = (wakeUpTime == nullTime ? sim.startTime : wakeUpTime)
	
	@assert(ambulance.index != nullIndex)
	@assert(ambulance.stationIndex != nullIndex)
	@assert(sim.startTime <= wakeUpTime)
	
	ambulance.status = ambSleeping
	# ambulance.stationIndex
	# ambulance.callIndex
	
	# create route that mimics ambulance driving from nowhere,
	# to a node (nearest to station), then to station, before simulation began
	ambulance.route = Route()
	ambulance.route.startLoc = Location()
	# ambulance.route.startTime = nullTime
	ambStation = sim.stations[ambulance.stationIndex]
	ambulance.route.endLoc = ambStation.location
	ambulance.route.endTime = sim.startTime
	ambulance.route.endFNode = ambStation.nearestNodeIndex
	
	# add wake up event
	addEvent!(sim.eventList; form = ambWakesUp, time = wakeUpTime, ambulance = ambulance)
end
