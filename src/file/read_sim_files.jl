# read simulation input files

function readAmbsFile(filename::String)
	tables = readTablesFromFile(filename)
	
	table = tables["ambulances"]
	n = size(table.data,1) # number of ambulances
	assert(n >= 1)
	
	# create ambulances from data in table
	ambulances = Vector{Ambulance}(n)
	columns = table.columns # shorthand
	for i = 1:n
		ambulances[i] = Ambulance()
		ambulances[i].index = columns["index"][i]
		ambulances[i].stationIndex = columns["stationIndex"][i]
		ambulances[i].class = AmbClass(columns["class"][i])
		assert(ambulances[i].index == i)
	end
	
	return ambulances
end

function readArcsFile(filename::String)
	tables = readTablesFromFile(filename)
	
	# get arcForm, numModes
	table = tables["miscData"]
	arcForm = table.columns["arcForm"][1] # "directed" or "undirected" arcs
	assert(in(arcForm, ["directed", "undirected"]))
	numModes = table.columns["numModes"][1]
	assert(numModes >= 1)
	
	# swap to arcs table
	table = tables["arcs"]
	
	# initialise arcs, travelTimes
	numArcs = size(table.data,1) # number of arcs, this may double if arcForm == "undirected"
	assert(numArcs >= 1)
	arcs = []
	travelTimes = []
	if arcForm == "directed"
		arcs = Vector{Arc}(numArcs)
		travelTimes = Array{Float,2}(numModes, numArcs)
	elseif arcForm == "undirected"
		arcs = Vector{Arc}(numArcs * 2)
		travelTimes = Array{Float,2}(numModes, numArcs * 2)
	end
	
	# check that table header has all travel times and modes
	for j = 1:numModes
		@assert(in(string("mode_", j), table.header), "Missing travel mode $j")
	end
	
	# fill arc fields with data in table
	c = table.columns # shorthand
	(indexCol = c["index"]); (fromNodeCol = c["fromNode"]); (toNodeCol = c["toNode"]) # shorthand, to avoid repeated dict lookups
	modeCols = [c[string("mode_", i)] for i = 1:numModes] # shorthand
	data = table.data # shorthand
	# get arc data from table
	for i = 1:numArcs
		arcs[i] = Arc()
		arcs[i].index = indexCol[i]
		arcs[i].fromNodeIndex = fromNodeCol[i]
		arcs[i].toNodeIndex = toNodeCol[i]
		# read travel times
		for j = 1:numModes
			travelTimes[j,i] = modeCols[j][i]
		end
		assert(arcs[i].index == i)
	end
	
	# if arcs are undirected, make them directed
	if arcForm == "undirected"
		# copy arcs, swap directions
		for i = 1:numArcs
			j = i + numArcs
			arcs[j] = Arc()
			arcs[j].index = j
			arcs[j].fromNodeIndex = arcs[i].toNodeIndex
			arcs[j].toNodeIndex = arcs[i].fromNodeIndex
		end
		# copy travel times
		travelTimes[:, (1:numArcs) + numArcs] = travelTimes[:, 1:numArcs]
	end
	
	# travel times should be positive
	assert(all(travelTimes .> 0))
	
	return arcs, travelTimes
end

function readCallsFile(filename::String)
	tables = readTablesFromFile(filename)
	
	# calls file has sim start time
	startTime = tables["miscData"].columns["startTime"][1]
	assert(startTime > nullTime)
	
	# create calls from data in calls table
	table = tables["calls"]
	data = table.data # shorthand
	n = size(data,1) # number of calls
	assert(n >= 1)
	c = table.columns # shorthand
	(indexCol = c["index"]); (priorityCol = c["priority"]); (xCol = c["x"]); (yCol = c["y"]); (arrivalTimeCol = c["arrivalTime"]); (dispatchDelayCol = c["dispatchDelay"]); (onSceneDurationCol = c["onSceneDuration"]); (transferDurationCol = c["transferDuration"]); (transferCol = c["transfer"]); (hospitalIndexCol = c["hospitalIndex"]) # shorthand, to avoid repeated dict lookups
	calls = Vector{Call}(n)
	for i = 1:n
		calls[i] = Call()
		calls[i].index = indexCol[i]
		calls[i].priority = Priority(priorityCol[i])
		calls[i].location.x = xCol[i]
		calls[i].location.y = yCol[i]
		calls[i].arrivalTime = arrivalTimeCol[i]
		calls[i].dispatchDelay = dispatchDelayCol[i]
		calls[i].onSceneDuration = onSceneDurationCol[i]
		calls[i].transferDuration = transferDurationCol[i]
		calls[i].transfer = transferCol[i]
		calls[i].hospitalIndex = hospitalIndexCol[i]
		
		assert(calls[i].index == i)
		assert(calls[i].priority != nullPriority)
		assert(calls[i].arrivalTime >= 0)
		assert(calls[i].dispatchDelay >= 0)
		assert(calls[i].onSceneDuration >= 0)
		if calls[i].transfer
			assert(calls[i].transferDuration >= 0)
		end
	end
	
	# check that calls are ordered by arrival time
	for i = 1:n-1
		assert(calls[i].arrivalTime <= calls[i+1].arrivalTime)
	end
	
	# check that first call is not before start time
	assert(startTime <= calls[1].arrivalTime)
	
	# check that either: no calls have same arrival time, or dispatch delay values are all > 0
	# This is to avoid potential bug where two calls arrive at same time with dispatch delay 0,
	# leading to calls being responded to in order of occurence in calls file, rather than by priority.
	assert(all(calls[i].arrivalTime < calls[i+1].arrivalTime for i = 1:n-1)
		|| all(calls[i].dispatchDelay > 0 for i = 1:n))
	
	return calls, startTime
end

function readCompTableFile(filename::String)
	tables = readTablesFromFile(filename)
	
	table = tables["compTable"]
	numAmbs = size(table.data,1)
	numStations = size(table.data,2) - 1
	columns = table.columns # shorthand
	
	# check that numAmbs column in table contains all ambulance counts
	assert(all(sort(columns["numAmbs"]) .== [1:numAmbs;]))
	
	# check that table header contains all stations
	for j = 1:numStations
		assert(in(string("station_", j), table.header))
	end
	
	# create compliance table
	compTable = Array{Int,2}(numAmbs, numStations)
	for i in columns["numAmbs"]
		for j = 1:numStations
			compTable[i,j] = columns[string("station_", j)][i]
		end
	end
	
	assert(all(compTable .>= 0))
	
	# check that compliance table rows add up to row indices
	assert(vec(sum(compTable, 2)) == [1:numAmbs;])
	
	return compTable
end

function readEventsFile(filename::String)
	tables = readTablesFromFile(filename)
	
	# input files
	table = tables["inputFiles"]
	inputFiles = table.columns["name"]
	fileChecksums = [parse(UInt, checksumString[2:end-1]) for checksumString in table.columns["checksum"]] # remove quote chars around number, convert to UInt
	
	# event dictionary
	table = tables["eventDict"]
	eventNames = table.columns["name"]
	eventKeys = table.columns["key"]
	n = length(eventKeys) # number of event types
	assert(length(eventNames) == n)
	assert(allunique(eventNames))
	assert(allunique(eventKeys))
	eventDict = Dict{Int,EventForm}()
	for i = 1:n
		eventDict[eventKeys[i]] = eval(parse(eventNames[i]))
	end
	
	# create events from data in events table
	table = tables["events"]
	c = table.columns # shorthand
	(eventKeyCol = c["eventKey"]); (timeCol = c["time"]); (ambIndexCol = c["ambIndex"]); (callIndexCol = c["callIndex"]); (stationIndexCol = c["stationIndex"]) # shorthand, to avoid repeated dict lookups
	data = table.data # shorthand
	fileEnded = (data[end,1] == "end") # true if writing all events to file ended before file was closed
	n = size(data,1) - fileEnded # number of events
	events = Vector{Event}(n)
	for i = 1:n
		events[i] = Event()
		events[i].form = eventDict[eventKeyCol[i]]
		events[i].time = timeCol[i]
		events[i].ambIndex = ambIndexCol[i]
		events[i].callIndex = callIndexCol[i]
		events[i].stationIndex = stationIndexCol[i]
		
		assert(events[i].form != nullEvent)
		assert(events[i].time > nullTime)
		assert(events[i].ambIndex != nullIndex || events[i].callIndex != nullIndex)
	end
	
	return events, fileEnded, inputFiles, fileChecksums
end

function readHospitalsFile(filename::String)
	tables = readTablesFromFile(filename)
	
	table = tables["hospitals"]
	n = size(table.data,1) # number of hospitals
	assert(n >= 1)
	
	# create hospitals from table data
	columns = table.columns # shorthand
	hospitals = Vector{Hospital}(n)
	for i = 1:n
		hospitals[i] = Hospital()
		hospitals[i].index = columns["index"][i]
		hospitals[i].location.x = columns["x"][i]
		hospitals[i].location.y = columns["y"][i]
		assert(hospitals[i].index == i)
	end
	
	return hospitals
end

function readMapFile(filename::String)
	tables = readTablesFromFile(filename)
	
	table = tables["map"]
	columns = table.columns # shorthand
	
	# create map from data in table
	map = Map()
	map.xMin = columns["xMin"][1]
	map.xMax = columns["xMax"][1]
	map.yMin = columns["yMin"][1]
	map.yMax = columns["yMax"][1]
	map.xScale = columns["xScale"][1]
	map.yScale = columns["yScale"][1]
	
	assert(map.xMin < map.xMax)
	assert(map.yMin < map.yMax)
	assert(0 < map.xScale)
	assert(0 < map.yScale)
	
	map.xRange = map.xMax - map.xMin
	map.yRange = map.yMax - map.yMin
	
	return map
end

function readNodesFile(filename::String)
	tables = readTablesFromFile(filename)
	
	table = tables["nodes"]
	n = size(table.data,1) # number of nodes
	assert(n >= 2)
	c = table.columns # shorthand
	(indexCol = c["index"]); (xCol = c["x"]); (yCol = c["y"]) # shorthand, to avoid repeated dict lookups
	accessCol = (haskey(c, "offRoadAccess") ? convert(Vector{Bool}, c["offRoadAccess"]) : fill(true, n))
	
	# create nodes from data in table
	data = table.data # shorthand
	nodes = Vector{Node}(n)
	for i = 1:n
		nodes[i] = Node()
		nodes[i].index = indexCol[i]
		nodes[i].location.x = xCol[i]
		nodes[i].location.y = yCol[i]
		nodes[i].offRoadAccess = accessCol[i]
		
		assert(nodes[i].index == i)
	end
	
	return nodes
end

function readPrioritiesFile(filename::String)
	tables = readTablesFromFile(filename)
	
	table = tables["priorities"]
	n = size(table.data,1) # number of priorities
	assert(n >= 1)
	assert(n <= 3) # have already hard-coded priorities: high, med, low; see defs.jl
	
	# read data from table
	columns = table.columns # shorthand
	targetResponseTimes = Vector{Float}(n)
	for i = 1:n
		assert(columns["priority"][i] == i)
		assert(eval(parse(columns["name"][i])) == Priority(i))
		targetResponseTimes[i] = columns["targetResponseTime"][i]
	end
	
	# target response times should be positive
	assert(all(targetResponseTimes .> 0))
	
	return targetResponseTimes
end

function readPriorityListFile(filename::String)
	tables = readTablesFromFile(filename)
	
	table = tables["priorityList"]
	n = size(table.data,1) # number of ambulances
	assert(n >= 1)
	
	# create priority list from data in table
	columns = table.columns # shorthand
	priorityList = Vector{Int}(n)
	for i = 1:n
		assert(columns["numAmbs"][i] == i)
		priorityList[i] = columns["stationIndex"][i]
	end
	
	# station indices should be positive
	assert(all(priorityList .>= 1))
	
	return priorityList
end

function readStationsFile(filename::String)
	tables = readTablesFromFile(filename)
	
	table = tables["stations"]
	n = size(table.data,1) # number of stations
	assert(n >= 1)
	
	# create stations from data in table
	columns = table.columns # shorthand
	stations = Vector{Station}(n)
	for i = 1:n
		stations[i] = Station()
		stations[i].index = columns["index"][i]
		stations[i].location.x = columns["x"][i]
		stations[i].location.y = columns["y"][i]
		stations[i].capacity = columns["capacity"][i]
		
		assert(stations[i].index == i)
		assert(stations[i].capacity >= 0)
	end
	
	return stations
end

function readTravelFile(filename::String)
	tables = readTablesFromFile(filename)
	
	# create travel object, will populate gradually
	travel = Travel()
	
	# travel modes
	table = tables["travelModes"]
	n = numModes = travel.numModes = size(table.data,1)
	assert(n >= 1)
	# create travel modes from data in table (network based data will need to be filled later)
	columns = table.columns # shorthand
	travel.modes = Vector{TravelMode}(n)
	for i = 1:n
		travel.modes[i] = TravelMode()
		travel.modes[i].index = columns["travelModeIndex"][i]
		travel.modes[i].offRoadSpeed = columns["offRoadSpeed"][i]
		
		assert(travel.modes[i].index == i)
		assert(travel.modes[i].offRoadSpeed > 0)
	end
	
	# travel sets
	table = tables["travelSets"]
	n = size(table.data,1)
	assert(n >= 1)
	# populate travel object from data in table
	columns = table.columns # shorthand
	numSets = travel.numSets = maximum(columns["travelSetIndex"])
	assert(sort(unique(columns["travelSetIndex"])) == [1:numSets;]) # all sets from 1:numSets should be used
	assert(sort(unique(columns["travelModeIndex"])) == [1:numModes;]) # all travel modes should be used
	travel.modeLookup = fill(nullIndex, numSets, 3) # number of priorities is fixed at 3, set defs.jl
	for i = 1:n
		setIndex = columns["travelSetIndex"][i]
		priority = eval(parse(columns["priority"][i]))
		modeIndex = columns["travelModeIndex"][i]
		
		assert(travel.modeLookup[setIndex, Int(priority)] == nullIndex) # value should not yet be filled
		assert(1 <= modeIndex <= numModes)
		travel.modeLookup[setIndex, Int(priority)] = modeIndex
	end
	assert(all(travel.modeLookup .!= nullIndex)) # check that travel.modeLookup is filled correctly
	
	# travel sets timing
	table = tables["travelSetsTiming"]
	# populate travel object from data in table
	columns = table.columns # shorthand
	travel.setsStartTimes = columns["startTime"]
	assert(travel.setsStartTimes[1] > nullTime)
	assert(issorted(travel.setsStartTimes, lt=<=)) # times should be strictly increasing
	travel.setsTimeOrder = columns["travelSetIndex"]
	assert(sort(unique(travel.setsTimeOrder)) == [1:numSets;]) # each travel set should be used at least once
	
	# misc
	travel.recentSetsStartTimesIndex = 1
	
	return travel
end

# read deployment policies from file
function readDeploymentPoliciesFile(filename::String)
	tables = readTablesFromFile(filename)
	
	# get counts
	table = tables["miscData"]
	numStations = table.columns["numStations"][1]
	numDepols = table.columns["numDepols"][1]
	
	# deployment policies
	table = tables["deploymentPolicies"]
	columns = table.columns # shorthand
	# check number of ambulances
	ambIndexCol = columns["ambIndex"]
	assert(ambIndexCol == collect(1:length(ambIndexCol)))
	# check that table header has all deployment policies
	for i = 1:numDepols
		@assert(in("policy $i, stationIndex", table.header), "Missing deployment policy $i")
	end
	# create deployment policies from table
	depols = Vector{Depol}(numDepols)
	for i = 1:numDepols
		depols[i] = columns["policy $i, stationIndex"]
		
		for value in depols[i]
			assert(in(value, 1:numStations))
		end
	end
	
	return depols, numStations
end
