# write sim object and output files

function writeAmbsFile(filename::String, ambulances::Vector{Ambulance})
	table = Table("ambulances", ["index", "stationIndex", "class"];
		rows = [[a.index, a.stationIndex, Int(a.class)] for a in ambulances])
	writeTablesToFile(filename, table)
end

function writeArcsFile(filename::String, arcs::Vector{Arc}, travelTimes::Array{Float}, arcForm::String)
	assert(arcForm == "directed" || arcForm == "undirected")
	numModes = size(travelTimes,2)
	miscTable = Table("miscData", ["arcForm", "numModes"]; rows = [[arcForm, numModes]])
	arcsTable = Table("arcs", vcat("index", "fromNode", "toNode", ["mode_$i" for i = 1:numModes]);
		rows = [vcat(a.index, a.fromNodeIndex, a.toNodeIndex, travelTimes[a.index,:]...) for a in arcs])
	writeTablesToFile(filename, [miscTable, arcsTable])
end

function writeCallsFile(filename::String, startTime::Float, calls::Vector{Call})
	miscTable = Table("miscData", ["startTime"]; rows = [[startTime]])
	callsTable = Table("calls", ["index", "priority", "x", "y", "arrivalTime", "dispatchDelay", "onSceneDuration", "transfer", "hospitalIndex", "transferDuration"];
		rows = [[c.index, Int(c.priority), c.location.x, c.location.y, c.arrivalTime, c.dispatchDelay, c.onSceneDuration, Int(c.transfer), c.hospitalIndex, c.transferDuration] for c in calls])
	writeTablesToFile(filename, [miscTable, callsTable])
end

function writeHospitalsFile(filename::String, hospitals::Vector{Hospital})
	table = Table("hospitals", ["index", "x", "y"];
		rows = [[h.index, h.location.x, h.location.y] for h in hospitals])
	writeTablesToFile(filename, table)
end

function writeMapFile(filename::String, map::Map)
	table = Table("map", ["xMin", "xMax", "yMin", "yMax", "xScale", "yScale"];
		rows = [[map.xMin, map.xMax, map.yMin, map.yMax, map.xScale, map.yScale]])
	writeTablesToFile(filename, table)
end

function writeNodesFile(filename::String, nodes::Vector{Node})
	table = Table("nodes", ["index", "x", "y", "offRoadAccess"];
		rows = [[n.index, n.location.x, n.location.y, Int(n.offRoadAccess)] for n in nodes])
	writeTablesToFile(filename, table)
end

function writePrioritiesFile(filename::String, targetResponseTimes::Vector{Float})
	table = Table("priorities", ["priority", "name", "targetResponseTime"];
		rows = [[i, string(Priority(i)), targetResponseTimes[i]] for i = 1:length(targetResponseTimes)])
	writeTablesToFile(filename, table)
end

function writeStationsFile(filename::String, stations::Vector{Station})
	table = Table("stations", ["index", "x", "y", "capacity"];
		rows = [[s.index, s.location.x, s.location.y, s.capacity] for s in stations])
	writeTablesToFile(filename, table)
end

function writeTravelFile(filename::String, travel::Travel)
	# travel modes table
	travelModesTable = Table("travelModes", ["travelModeIndex", "offRoadSpeed"];
		rows = [[tm.index, tm.offRoadSpeed] for tm in travel.modes])
	
	# travel sets table
	tml = travel.modeLookup # shorthand
	assert(size(tml) == (length(travel.modes), 3)) # should have value for each combination of travel mode and priority
	travelSetsTable = Table("travelSets", ["travelSetIndex", "priority", "travelModeIndex"];
		rows = [[ind2sub(tml,i)[1], string(Priority(ind2sub(tml,i)[2])), tml[i]] for i = 1:length(tml)])
	
	# travel sets timing table
	startTimes = travel.setsStartTimes # shorthand
	setIndices = travel.setsTimeOrder # shorthand
	travelSetsTimingTable = Table("travelSetsTiming", ["startTime", "travelSetIndex"];
		rows = [[startTimes[i], setIndices[i]] for i = 1:length(startTimes)])
	
	writeTablesToFile(filename, [travelModesTable, travelSetsTable, travelSetsTimingTable])
end

# opens output files for writing during simulation
# note: should have field sim.resim.use set/fixed before calling this function
function openOutputFiles!(sim::Simulation)
	if !sim.writeOutput; return; end
	
	println("output path: ", sim.outputPath)
	outputFilePath(name::String) = sim.outputFiles[name].path
	
	# create output path if it does not already exist
	if !isdir(sim.outputPath)
		mkdir(sim.outputPath)
	end
	
	# open output files for writing
	# currently, only need to open events file
	if !sim.resim.use # otherwise, existing events file is used for resimulation
		sim.outputFiles["events"].iostream = open(outputFilePath("events"), "w")
		sim.eventsFileIO = sim.outputFiles["events"].iostream # shorthand
		
		# save checksum of input files
		inputFiles = sort([name for (name, file) in sim.inputFiles])
		fileChecksumStrings = [string("'", sim.inputFiles[name].checksum, "'") for name in inputFiles]
		writeTablesToFile!(sim.eventsFileIO, Table("inputFiles", ["name", "checksum"]; cols = [inputFiles, fileChecksumStrings]))
		
		# save events with a key, to reduce file size
		eventForms = instances(EventForm)
		eventKeys = [Int(eventForm) for eventForm in eventForms]
		eventNames = [string(eventForm) for eventForm in eventForms]
		writeTablesToFile!(sim.eventsFileIO, Table("eventDict", ["key", "name"]; cols = [eventKeys, eventNames]))
		
		# write events table name and header
		writeDlmLine!(sim.eventsFileIO, "events")
		writeDlmLine!(sim.eventsFileIO, "time", "eventKey", "ambIndex", "callIndex", "stationIndex")
	end
end

function closeOutputFiles!(sim::Simulation)
	if !sim.writeOutput; return; end
	
	if !sim.resim.use
		writeDlmLine!(sim.eventsFileIO, "end")
		close(sim.eventsFileIO)
	end
end

function writeEventToFile!(sim::Simulation, event::Event)
	if !sim.writeOutput || sim.resim.use; return; end
	
	# find station index of ambulance (if any)
	stationIndex = nullIndex
	if event.ambIndex != nullIndex
		stationIndex = sim.ambulances[event.ambIndex].stationIndex
	end
	
	writeDlmLine!(sim.eventsFileIO, @sprintf("%0.5f", event.time), Int(event.form),
		event.ambIndex, event.callIndex, stationIndex)
		
	# flush(sim.eventsFileIO)
end

function writeStatsFiles!(sim::Simulation)
	if !sim.writeOutput; return; end
	
	timeRounding = 6 # number of digits to keep in time values, after decimal place
	
	outputFilePath(name::String) = sim.outputFiles[name].path
	
	println("saving stats to folder: ", sim.outputPath)
	
	# save ambulance stats
	writeTablesToFile(outputFilePath("ambulances"), Table("ambStats",
		["index", "stationIndex", "totalTravelTime", "totalBusyTime", "numCallsTreated", "numCallsTransferred", "atStationDispatches", "onRoadDispatches", "afterServiceDispatches", "numDiversions"];
		rows = [vcat(a.index, a.stationIndex, round.([a.totalTravelTime, a.totalBusyTime], timeRounding)..., a.numCallsTreated, a.numCallsTransferred, a.atStationDispatches, a.onRoadDispatches, a.afterServiceDispatches, a.numDiversions) for a in sim.ambulances]))
	
	# save call stats
	writeTablesToFile(outputFilePath("calls"), Table("callStats",
		["index", "priority", "ambIndex", "transfer", "arrivalTime", "dispatchDelay", "onSceneDuration", "hospitalArrivalTime", "transferDuration", "dispatchTime", "ambArrivalTime", "hospitalIndex", "numBumps", "wasQueued"];
		rows = [vcat(c.index, Int(c.priority), c.ambIndex, Int(c.transfer), round.([c.arrivalTime, c.dispatchDelay, c.onSceneDuration, c.hospitalArrivalTime, c.transferDuration * c.transfer, c.dispatchTime, c.ambArrivalTime], timeRounding)..., c.hospitalIndex, c.numBumps, Int(c.wasQueued)) for c in sim.calls]))
	
	# save hospital stats file
	writeTablesToFile(outputFilePath("hospitals"), Table("hospitalStats", ["index", "numTransfers"];
		rows = [[h.index, h.numTransfers] for h in sim.hospitals]))
end

# write deployment policies to file
function writeDeploymentPoliciesFile(filename::String, depols::Vector{Depol}, numStations::Int)
	numAmbs = length(depols[1])
	assert(numStations >= maximum([maximum(d) for d in depols]))
	numDepols = length(depols)
	
	miscTable = Table("miscData", ["numStations", "numDepols"]; rows = [[numStations, numDepols]])
	deploymentPoliciesTable = Table("deploymentPolicies",
		["ambIndex", ["policy $i, stationIndex" for i = 1:numDepols]...];
		cols = [collect(1:numAmbs), depols...])
	writeTablesToFile(filename, [miscTable, deploymentPoliciesTable])
end
