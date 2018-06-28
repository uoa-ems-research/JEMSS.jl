# run configuration xml file
function runConfig(configFilename::String)
	# read input files, initialise simulation
	sim = initSimulation(configFilename; allowResim = true, createBackup = false, allowWriteOutput = true)
	
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
function initSimulation(configFilename::String;
	allowResim::Bool = false, createBackup::Bool = true, allowWriteOutput::Bool = false)
	
	# read sim config xml file
	rootElt = xmlFileRoot(configFilename)
	@assert(name(rootElt) == "simConfig", string("xml root has incorrect name: ", name(rootElt)))
	
	# for progress messages:
	t = Vector{Float}(1)
	initMessage(t, msg) = (t[1] = time(); print(msg))
	initTime(t) = println(": ", round(time() - t[1], 2), " seconds")
	
	##################
	# sim config
	
	initMessage(t, "reading config file data")
	
	sim = Simulation()
	sim.configRootElt = rootElt
	
	# input
	sim.inputPath = abspath(eltContentInterpVal(rootElt, "inputPath"))
	simFilesElt = findElt(rootElt, "simFiles")
	inputFiles = childrenNodeNames(simFilesElt)
	sim.inputFiles = Dict{String,File}()
	for inputFile in inputFiles
		file = File()
		file.name = eltContent(simFilesElt, inputFile)
		file.path = joinpath(sim.inputPath, file.name)
		if inputFile != "rNetTravels" # do not need checksum of rNetTravels file
			file.checksum = fileChecksum(file.path)
		end
		sim.inputFiles[inputFile] = file
	end
	
	# output
	sim.writeOutput = allowWriteOutput && eltContentVal(rootElt, "writeOutput")
	sim.outputPath = abspath(eltContentInterpVal(rootElt, "outputPath"))
	outputFilesElt = findElt(rootElt, "outputFiles")
	outputFiles = childrenNodeNames(outputFilesElt)
	sim.outputFiles = Dict{String,File}()
	for outputFile in outputFiles
		file = File()
		file.name = eltContent(outputFilesElt, outputFile)
		file.path = joinpath(sim.outputPath, file.name)
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
	sim.targetResponseTimes = readPrioritiesFile(simFilePath("priorities"))
	sim.travel = readTravelFile(simFilePath("travel"))
	
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
	println("fNodes: ", length(net.fGraph.nodes), ", rNodes: ", length(net.rGraph.nodes))
	
	initMessage(t, "checking rGraph")
	checkGraph(net.rGraph, map)
	initTime(t)
	
	if rNetTravelsLoaded != []
		println("using data from rNetTravels file")
		try
			initMessage(t, "creating rNetTravels from fNetTravels")
			createRNetTravelsFromFNetTravels!(net; rNetTravelsLoaded = rNetTravelsLoaded)
			initTime(t)
		catch
			println()
			warn("failed to use data from rNetTravels file")
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
	assert(travel.setsStartTimes[1] <= sim.startTime)
	assert(length(net.fNetTravels) == travel.numModes)
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
	
	println("nodes: ", length(fGraph.nodes), ", grid size: ", nx, " x ", ny)
	
	##################
	# sim - ambulances, calls, hospitals, stations...
	
	initMessage(t, "adding ambulances, calls, etc")
	
	# for each call, hospital, and station, find neareset node
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
	sim.eventList = Vector{Event}(0)
	
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
		fNetTravel.fNodeNearestHospitalIndex = Vector{Int}(numFNodes)
		travelModeIndex = fNetTravel.modeIndex # shorthand
		travelMode = travel.modes[travelModeIndex] # shorthand
		for node in fGraph.nodes
			# find nearest hospital to node
			minTime = Inf
			nearestHospitalIndex = nullIndex
			for hospital in sim.hospitals
				(travelTime, rNodes) = shortestPathTravelTime(net, travelModeIndex, node.index, hospital.nearestNodeIndex)
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
	
	# move up
	mud = sim.moveUpData # shorthand
	moveUpElt = findElt(decisionElt, "moveUp")
	moveUpModuleName = eltContent(moveUpElt, "module")
	
	if moveUpModuleName == "none"
		mud.useMoveUp = false
		mud.moveUpModule = nullMoveUpModule
		println("not using move up")
	else
		initMessage(t, "initialising move up")
		
		mud.useMoveUp = true
		
		if moveUpModuleName == "comp_table"
			mud.moveUpModule = compTableModule
			compTableElt = findElt(moveUpElt, "compTable")
			compTableFilename = eltContent(compTableElt, "filename")
			initCompTable!(sim, joinpath(sim.inputPath, compTableFilename))
			
		elseif moveUpModuleName == "dmexclp"
			mud.moveUpModule = dmexclpModule
			dmexclpElt = findElt(moveUpElt, "dmexclp")
			initDmexclp!(sim;
				coverTime = eltContentVal(dmexclpElt, "coverTime"),
				coverTravelPriority = eltContentVal(dmexclpElt, "coverTravelPriority"),
				busyFraction = eltContentVal(dmexclpElt, "busyFraction"),
				demandRasterFilename = eltContent(dmexclpElt, "demandRasterFilename"))
			
		elseif moveUpModuleName == "priority_list"
			mud.moveUpModule = priorityListModule
			priorityListElt = findElt(moveUpElt, "priorityList")
			priorityListFilename = eltContent(priorityListElt, "filename")
			initPriorityList!(sim, joinpath(sim.inputPath, priorityListFilename))
			
		elseif moveUpModuleName == "zhang_ip"
			mud.moveUpModule = zhangIpModule
			zhangIpElt = findElt(moveUpElt, "zhangIp")
			initZhangIp!(sim;
				busyFraction = eltContentVal(zhangIpElt, "busyFraction"),
				travelTimeCost = eltContentVal(zhangIpElt, "travelTimeCost"),
				maxIdleAmbTravelTime = eltContentVal(zhangIpElt, "maxIdleAmbTravelTime"),
				maxNumNearestStations = eltContentVal(zhangIpElt, "maxNumNearestStations"))
			
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
		
		println("using move up module: ", moveUpModuleName)
	end
	
	##################
	# resimulation
	
	if allowResim
		sim.resim.use = eltContentVal(rootElt, "resim")
		if sim.resim.use
			initMessage(t, "")
			initResimulation!(sim)
			print("initialised resimulation")
			initTime(t)
		end
	end
	
	##################
	# sim backup
	
	if createBackup
		initMessage(t, "creating sim backup")
		backupSim!(sim) # for restarting sim
		initTime(t)
	end
	
	return sim
end

# initialise given ambulance
# sets ambulance as sleeping, creates wake up event
function initAmbulance!(sim::Simulation, ambulance::Ambulance;
	wakeUpTime::Float = nullTime)
	wakeUpTime = (wakeUpTime == nullTime ? sim.startTime : wakeUpTime)
	
	assert(ambulance.index != nullIndex)
	assert(ambulance.stationIndex != nullIndex)
	assert(sim.startTime <= wakeUpTime)
	
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
