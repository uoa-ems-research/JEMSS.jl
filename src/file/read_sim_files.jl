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

# read simulation input files

function readAmbsFile(filename::String)
	tables = readTablesFromFile(filename)
	
	table = tables["ambulances"]
	n = size(table.data,1) # number of ambulances
	@assert(n >= 1)
	
	# create ambulances from data in table
	ambulances = Vector{Ambulance}(undef, n)
	columns = table.columns # shorthand
	for i = 1:n
		ambulances[i] = Ambulance()
		ambulances[i].index = columns["index"][i]
		ambulances[i].stationIndex = columns["stationIndex"][i]
		ambulances[i].class = AmbClass(columns["class"][i])
		@assert(ambulances[i].index == i)
	end
	
	return ambulances
end

function readArcsFile(filename::String;
	keepAllFields::Bool = false)
	tables = readTablesFromFile(filename)
	
	# get arcForm, numModes
	table = tables["miscData"]
	arcForm = table.columns["arcForm"][1] # "directed" or "undirected" arcs
	@assert(in(arcForm, ["directed", "undirected"]))
	numModes = table.columns["numModes"][1]
	@assert(numModes >= 1)
	
	# swap to arcs table
	table = tables["arcs"]
	
	# initialise arcs, travelTimes
	numArcs = size(table.data,1) # number of arcs, this may double if arcForm == "undirected"
	@assert(numArcs >= 1)
	arcs = []
	travelTimes = []
	if arcForm == "directed"
		arcs = Vector{Arc}(undef, numArcs)
		travelTimes = Array{Float,2}(undef, numModes, numArcs)
	elseif arcForm == "undirected"
		arcs = Vector{Arc}(undef, numArcs * 2)
		travelTimes = Array{Float,2}(undef, numModes, numArcs * 2)
	end
	
	# check that table header has all travel times and modes
	for j = 1:numModes
		@assert(in(string("mode_", j), table.header), "Missing travel mode $j")
	end
	
	# fill arc fields with data in table
	c = table.columns # shorthand
	(indexCol = c["index"]); (fromNodeCol = c["fromNode"]); (toNodeCol = c["toNode"]) # shorthand, to avoid repeated dict lookups
	(hasDistCol, distCol) = haskey(c, "distance") ? (true, c["distance"]) : (false, fill(NaN, numArcs)) # shorthand
	modeCols = [c[string("mode_", i)] for i = 1:numModes] # shorthand
	fieldNames = setdiff(table.header, keepAllFields ? [] : ["index", "fromNode", "toNode", "distance", [string("mode_", i) for i = 1:numModes]...])
	rowsFields = tableRowsFieldDicts(table, fieldNames)
	# get arc data from table
	for i = 1:numArcs
		arcs[i] = Arc()
		arcs[i].index = indexCol[i]
		arcs[i].fromNodeIndex = fromNodeCol[i]
		arcs[i].toNodeIndex = toNodeCol[i]
		arcs[i].distance = distCol[i]
		arcs[i].fields = rowsFields[i]
		# read travel times
		for j = 1:numModes
			travelTimes[j,i] = modeCols[j][i]
		end
		@assert(arcs[i].index == i)
		if hasDistCol @assert(arcs[i].distance >= 0) end
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
			arcs[j].distance = arcs[i].distance
			arcs[j].fields = deepcopy(arcs[i].fields)
		end
		# copy travel times
		travelTimes[:, (1:numArcs) .+ numArcs] = travelTimes[:, 1:numArcs]
	end
	
	# travel times should be positive
	@assert(all(travelTimes .> 0))
	
	return arcs, travelTimes
end

function readCallsFile(filename::String)
	tables = readTablesFromFile(filename)
	
	# calls file has sim start time
	startTime = tables["miscData"].columns["startTime"][1]
	@assert(startTime > nullTime)
	
	# create calls from data in calls table
	table = tables["calls"]
	data = table.data # shorthand
	n = size(data,1) # number of calls
	@assert(n >= 1)
	c = table.columns # shorthand
	(indexCol = c["index"]); (priorityCol = c["priority"]); (xCol = c["x"]); (yCol = c["y"]); (arrivalTimeCol = c["arrivalTime"]); (dispatchDelayCol = c["dispatchDelay"]); (onSceneDurationCol = c["onSceneDuration"]); (hospitalIndexCol = c["hospitalIndex"]) # shorthand, to avoid repeated dict lookups
	if haskey(c, "transport") && haskey(c, "handoverDuration")
		(transportCol = c["transport"]); (handoverDurationCol = c["handoverDuration"]);
	else # compat
		@assert(haskey(c, "transfer") && haskey(c, "transferDuration"))
		@warn("`transfer` and `transferDuration` headers in calls file are deprecated, use `transport` and `handoverDuration` instead.")
		(transportCol = c["transfer"]); (handoverDurationCol = c["transferDuration"]);
	end
	
	calls = Vector{Call}(undef, n)
	for i = 1:n
		calls[i] = Call()
		calls[i].index = indexCol[i]
		calls[i].priority = Priority(Int(priorityCol[i]))
		calls[i].location.x = xCol[i]
		calls[i].location.y = yCol[i]
		calls[i].arrivalTime = arrivalTimeCol[i]
		calls[i].dispatchDelay = dispatchDelayCol[i]
		calls[i].onSceneDuration = onSceneDurationCol[i]
		calls[i].transport = transportCol[i]
		calls[i].handoverDuration = calls[i].transport ? handoverDurationCol[i] : 0.0
		calls[i].hospitalIndex = hospitalIndexCol[i]
		
		@assert(calls[i].index == i)
		@assert(calls[i].priority != nullPriority)
		@assert(calls[i].arrivalTime >= 0)
		@assert(calls[i].dispatchDelay >= 0)
		@assert(calls[i].onSceneDuration >= 0)
		@assert(calls[i].handoverDuration >= 0)
	end
	
	# check that calls are ordered by arrival time
	for i = 1:n-1
		@assert(calls[i].arrivalTime <= calls[i+1].arrivalTime)
	end
	
	# check that first call is not before start time
	@assert(startTime <= calls[1].arrivalTime)
	
	# check that either: no calls have same arrival time, or dispatch delay values are all > 0
	# This is to avoid potential bug where two calls arrive at same time with dispatch delay 0,
	# leading to calls being responded to in order of occurence in calls file, rather than by priority.
	@assert(all(calls[i].arrivalTime < calls[i+1].arrivalTime for i = 1:n-1)
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
	@assert(all(sort(columns["numAmbs"]) .== [1:numAmbs;]))
	
	# check that table header contains all stations
	for j = 1:numStations
		@assert(in(string("station_", j), table.header))
	end
	
	# create compliance table
	compTable = CompTable(undef,numAmbs, numStations)
	for i in columns["numAmbs"]
		for j = 1:numStations
			compTable[i,j] = columns[string("station_", j)][i]
		end
	end
	
	@assert(all(compTable .>= 0))
	
	# check that compliance table rows add up to row indices
	@assert(vec(sum(compTable, dims=2)) == [1:numAmbs;])
	
	return compTable
end

function readDemandFile(filename::String)
	tables = readTablesFromFile(filename)
	
	# create demand object, will populate gradually
	demand = Demand()
	
	# demand rasters
	table = tables["demandRasters"]
	n = numRasters = demand.numRasters = size(table.data,1)
	@assert(n >= 1)
	columns = table.columns # shorthand
	demand.rasters = Vector{Raster}(undef, n)
	demand.rasterFilenames = Vector{String}(undef, n)
	allunique(columns["rasterFilename"]) || @warn("There are duplicate raster filenames in the demand file, this is unnecessary.")
	demandFileDir = dirname(realpath(filename))
	for i = 1:n
		@assert(columns["rasterIndex"][i] == i)
		
		rasterFilename = String(columns["rasterFilename"][i])
		rasterFilename = joinPathIfNotAbs(demandFileDir, interpolateString(rasterFilename)) # raster filename should be absolute (after interpolating), or relative to the demand file directory
		@assert(isfile(rasterFilename))
		demand.rasterFilenames[i] = rasterFilename
		demand.rasters[i] = readRasterFile(rasterFilename)
	end
	
	# demand modes
	table = tables["demandModes"]
	n = numModes = demand.numModes = size(table.data,1)
	@assert(n >= 1)
	columns = table.columns # shorthand
	demand.modes = Vector{DemandMode}(undef, n)
	for i = 1:n
		@assert(columns["modeIndex"][i] == i)
		
		demand.modes[i] = DemandMode()
		demand.modes[i].index = i
		rasterIndex = demand.modes[i].rasterIndex = columns["rasterIndex"][i]
		@assert(1 <= rasterIndex <= numRasters)
		demand.modes[i].raster = demand.rasters[rasterIndex]
		demand.modes[i].priority = columns["priority"][i] |> Meta.parse |> eval
		demand.modes[i].arrivalRate = columns["arrivalRate"][i]
		@assert(demand.modes[i].arrivalRate >= 0)
		demand.modes[i].rasterMultiplier = demand.modes[i].arrivalRate / sum(demand.modes[i].raster.z)
	end
	
	# demand sets
	table = tables["demandSets"]
	n = size(table.data,1)
	@assert(n >= 1)
	columns = table.columns # shorthand
	numSets = demand.numSets = maximum(columns["setIndex"])
	@assert(sort(unique(columns["setIndex"])) == [1:numSets;]) # all sets from 1:numSets should be used
	demand.modeLookup = fill(nullIndex, numSets, numPriorities)
	for i = 1:n
		setIndex = columns["setIndex"][i]
		@assert(setIndex == i)
		modeIndices = columns["modeIndices"][i] |> Meta.parse |> eval
		for modeIndex in modeIndices
			@assert(1 <= modeIndex <= numModes)
			priorityIndex = Int(demand.modes[modeIndex].priority)
			@assert(demand.modeLookup[setIndex,priorityIndex] == nullIndex) # value should not yet be filled
			demand.modeLookup[setIndex,priorityIndex] = modeIndex
		end
	end
	@assert(all(demand.modeLookup .!= nullIndex)) # check that demand.modeLookup is filled correctly
	# have not checked that all modes from 1:numModes are used
	
	# demand sets timing
	table = tables["demandSetsTiming"]
	columns = table.columns # shorthand
	demand.setsStartTimes = columns["startTime"]
	@assert(demand.setsStartTimes[1] > nullTime)
	@assert(issorted(demand.setsStartTimes, lt=<=)) # times should be strictly increasing
	demand.setsTimeOrder = columns["setIndex"]
	@assert(sort(unique(demand.setsTimeOrder)) == [1:numSets;]) # each demand set should be used at least once
	
	# misc
	demand.recentSetsStartTimesIndex = 1
	
	return demand
end

function readDemandCoverageFile(filename::String)
	tables = readTablesFromFile(filename)
	
	dc = demandCoverage = DemandCoverage()
	
	# demand cover times
	table = tables["coverTimes"]
	n = size(table.data,1) # number of priorities
	@assert(n == numPriorities)
	columns = table.columns # shorthand
	for i = 1:n
		priority = eval(Meta.parse(columns["demandPriority"][i]))
		coverTime = dc.coverTimes[priority] = columns["coverTime"][i]
		@assert(coverTime >= 0)
	end
	@assert(all(p -> haskey(dc.coverTimes, p), priorities), "coverTimes not set for all priorities")
	
	# number of points to represent demand per raster cell
	table = tables["demandRasterCellNumPoints"]
	@assert(size(table.data,1) == 1) # should be single row
	columns = table.columns # shorthand
	rows = dc.rasterCellNumRows = columns["rows"][1]
	cols = dc.rasterCellNumCols = columns["cols"][1]
	@assert(rows >= 1)
	@assert(cols >= 1)
	
	return dc
end

function readEventsFile(filename::String)
	tables = readTablesFromFile(filename)
	
	# input files
	table = tables["inputFiles"]
	inputFiles = table.columns["name"]
	fileChecksums = [parse(UInt, checksumString[2:end-1]) for checksumString in table.columns["checksum"]] # remove quote chars around number, convert to UInt
	
	# event dictionary and filter
	table = tables["eventDict"]
	eventNames = table.columns["name"]
	eventKeys = table.columns["key"]
	n = length(eventKeys) # number of event types
	filter = haskey(table.columns, "filter") ? Bool.(table.columns["filter"]) : fill(true,n)
	@assert(length(eventNames) == n)
	@assert(allunique(eventNames))
	@assert(allunique(eventKeys))
	eventDict = Dict{Int,EventForm}()
	eventFilter = Dict([e => true for e in instances(EventForm)]) # eventFilter[eventForm] = true if events of eventForm are included in file, false otherwise
	for i = 1:n
		eventForm = eval(Meta.parse(eventNames[i]))
		eventDict[eventKeys[i]] = eventForm
		eventFilter[eventForm] = filter[i]
	end
	
	# create events from data in events table
	table = tables["events"]
	c = table.columns # shorthand
	(eventIndexCol = c["index"]); (parentIndexCol = c["parentIndex"]); (eventKeyCol = c["eventKey"]); (timeCol = c["time"]); (ambIndexCol = c["ambIndex"]); (callIndexCol = c["callIndex"]); (stationIndexCol = c["stationIndex"]) # shorthand, to avoid repeated dict lookups
	data = table.data # shorthand
	fileEnded = (data[end,1] == "end") # true if writing all events to file ended before file was closed
	n = size(data,1) - fileEnded # number of events in events table, can be less than actual number of events simulated
	events = Vector{Event}(undef, n)
	for i = 1:n
		event = events[i] = Event()
		event.index = eventIndexCol[i]
		event.parentIndex = parentIndexCol[i]
		event.form = eventDict[eventKeyCol[i]]
		event.time = timeCol[i]
		event.ambIndex = ambIndexCol[i]
		event.callIndex = callIndexCol[i]
		event.stationIndex = stationIndexCol[i]
		
		@assert(event.index > event.parentIndex)
		if i > 1 @assert(event.index > events[i-1].index) end
		@assert(event.form != nullEvent)
		@assert(eventFilter[event.form])
		@assert(event.time > nullTime)
		@assert(event.ambIndex != nullIndex || event.callIndex != nullIndex)
	end
	
	# set eventsChildren; eventsChildren[i] gives events that are children of the ith event that occurred
	eventsChildren = [Vector{Event}() for i = 1:events[end].index]
	for event in events
		j = event.parentIndex
		if j != nullIndex
			push!(eventsChildren[j], event)
		end
	end
	
	return events, eventsChildren, eventFilter, fileEnded, inputFiles, fileChecksums
end

# read geographic data from file and apply f to data
function readGeoFile(f::Function, filename::String)
	if pkgVersions["ArchGDAL"] >= v"0.3"
		return ArchGDAL.read(filename) do dataset
			f(dataset)
		end
	else
		return ArchGDAL.registerdrivers() do
			ArchGDAL.read(filename) do dataset
				f(dataset)
			end
		end
	end
end

function readHospitalsFile(filename::String)
	tables = readTablesFromFile(filename)
	
	table = tables["hospitals"]
	n = size(table.data,1) # number of hospitals
	@assert(n >= 1)
	
	# create hospitals from table data
	columns = table.columns # shorthand
	hospitals = Vector{Hospital}(undef, n)
	for i = 1:n
		hospitals[i] = Hospital()
		hospitals[i].index = columns["index"][i]
		hospitals[i].location.x = columns["x"][i]
		hospitals[i].location.y = columns["y"][i]
		@assert(hospitals[i].index == i)
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
	
	@assert(map.xMin < map.xMax)
	@assert(map.yMin < map.yMax)
	@assert(0 < map.xScale)
	@assert(0 < map.yScale)
	
	map.xRange = map.xMax - map.xMin
	map.yRange = map.yMax - map.yMin
	
	return map
end

function readNodesFile(filename::String)
	tables = readTablesFromFile(filename)
	
	table = tables["nodes"]
	n = size(table.data,1) # number of nodes
	@assert(n >= 2)
	c = table.columns # shorthand
	(indexCol = c["index"]); (xCol = c["x"]); (yCol = c["y"]) # shorthand, to avoid repeated dict lookups
	accessCol = (haskey(c, "offRoadAccess") ? convert(Vector{Bool}, c["offRoadAccess"]) : fill(true, n))
	fieldNames = setdiff(table.header, ["index", "x", "y", "offRoadAccess"])
	rowsFields = tableRowsFieldDicts(table, fieldNames)
	
	# create nodes from data in table
	nodes = Vector{Node}(undef, n)
	for i = 1:n
		nodes[i] = Node()
		nodes[i].index = indexCol[i]
		nodes[i].location.x = xCol[i]
		nodes[i].location.y = yCol[i]
		nodes[i].offRoadAccess = accessCol[i]
		nodes[i].fields = rowsFields[i]
		
		@assert(nodes[i].index == i)
	end
	
	return nodes
end

function readPrioritiesFile(filename::String)
	tables = readTablesFromFile(filename)
	
	table = tables["priorities"]
	n = size(table.data,1) # number of priorities
	@assert(n >= 1)
	@assert(n <= numPriorities)
	
	# read data from table
	columns = table.columns # shorthand
	targetResponseDurationString = haskey(columns, "targetResponseDuration") ? "targetResponseDuration" : "targetResponseTime" # compat for "targetResponseTime"
	targetResponseDurations = Vector{Float}(undef, n)
	responseTravelPriorities = Dict([p => p for p in priorities]) # default is to have response travel priority equal to call priority
	for i = 1:n
		@assert(columns["priority"][i] == i)
		@assert(eval(Meta.parse(columns["name"][i])) == Priority(i))
		targetResponseDurations[i] = columns[targetResponseDurationString][i]
		
		if haskey(columns, "responseTravelPriority")
			responseTravelPriorities[Priority(i)] = eval(Meta.parse(columns["responseTravelPriority"][i]))
		end
	end
	
	# target response durations should be positive
	@assert(all(targetResponseDurations .> 0))
	
	return targetResponseDurations, responseTravelPriorities
end

function readPriorityListFile(filename::String)::PriorityList
	tables = readTablesFromFile(filename)
	
	table = tables["priorityList"]
	n = size(table.data,1) # number of ambulances
	@assert(n >= 1)
	
	# create priority list from data in table
	columns = table.columns # shorthand
	itemCol = haskey(columns, "item") ? columns["item"] : columns["numAmbs"] # allow old column name "numAmbs"
	priorityList = PriorityList(undef,n)
	for i = 1:n
		@assert(itemCol[i] == i)
		priorityList[i] = columns["stationIndex"][i]
	end
	
	# station indices should be positive
	@assert(all(priorityList .>= 1))
	
	return priorityList
end

function readPriorityListsFile(filename::String)::Vector{PriorityList}
	tables = readTablesFromFile(filename)
	
	table = tables["priorityLists"]
	(numAmbs, n) = size(table.data)
	@assert(numAmbs >= 1)
	numPriorityLists = n - 1
	@assert(numPriorityLists >= 1) # number of priority lists
	
	# create priority list from data in table
	columns = table.columns # shorthand
	itemCol = haskey(columns, "item") ? columns["item"] : columns["numAmbs"] # allow old column name "numAmbs"
	priorityLists = PriorityList[]
	for i = 1:numPriorityLists
		priorityList = PriorityList(undef,numAmbs)
		for j = 1:numAmbs
			@assert(itemCol[j] == j)
			priorityList[j] = columns["priorityList_$i stationIndex"][j]
		end
		
		# station indices should be positive
		@assert(all(priorityList .>= 1))
		
		push!(priorityLists, priorityList)
	end
	
	return priorityLists
end


# read raster file using ArchGDAL package, return as custom Raster type
function readRasterFile(rasterFilename::String)
	@assert(isfile(rasterFilename))
	
	# open raster, get data, close raster
	# only gets first raster band, may change this later
	(geoTransform, z) = readGeoFile(rasterFilename) do dataset
		rasterband = ArchGDAL.getband(dataset, 1)
		return (ArchGDAL.getgeotransform(dataset), ArchGDAL.read(rasterband))
	end
	
	# shorthand:
	(nx, ny) = size(z) # number of cells in x and y directions
	(x1, dxdi, dxdj, y1, dydi, dydj) = geoTransform # dxdi is change in x per change in index i, for z[i,j]; likewise for dxdj, dydi, dydj
	dx = dxdi # width of cells in x direction
	dy = dydj # height of cells in y direction (may be negative)
	
	# data checks
	@assert(dxdj == 0 && dydi == 0) # otherwise raster is sloping, so changing x index changes y value, and vice-versa
	@assert(dx > 0)
	@assert(dy != 0)
	
	# convert data for easier use
	# find x and y vectors to represent raster cell centres, make sure values are increasing
	xMin = x1 + 0.5*dx # we know dx > 0
	# (xMin, dx, z) = (dx > 0) ? (x1 + 0.5*dx, dx, z) : (x1 + (nx-0.5)*dx, -dx, reverse(z, dims = 1))
	(yMin, dy, z) = (dy > 0) ? (y1 + 0.5*dy, dy, z) : (y1 + (ny-0.5)*dy, -dy, reverse(z, dims = 2))
	x = collect(range(xMin, step=dx, length=nx))
	y = collect(range(yMin, step=dy, length=ny))
	
	return Raster(x, y, z)
end

function readRedispatchFile(filename::String)::Redispatch
	tables = readTablesFromFile(filename)
	
	redispatch = Redispatch()
	
	table = tables["miscData"]
	redispatch.allow = Bool(table.columns["allowRedispatch"][1])
	
	table = tables["redispatchConditions"]
	n = size(table.data,1)
	priorityPairUsed = fill(false, numPriorities, numPriorities)
	for i = 1:n
		p1 = eval(Meta.parse(table.columns["fromCallPriority"][i]))
		p2 = eval(Meta.parse(table.columns["toCallPriority"][i]))
		redispatch.conditions[Int(p1),Int(p2)] = Bool(table.columns["allowRedispatch"][i])
		
		# check that condition will not overwrite existing condition
		@assert(priorityPairUsed[Int(p1),Int(p2)] == false, "Redispatch file has multiple conditions for call priority pair [$p1, $p2].")
		priorityPairUsed[Int(p1),Int(p2)] = true
		
		p1 == p2 && @assert(!redispatch.conditions[Int(p1),Int(p2)], "Redispatch file has condition that allows redispatch for calls of same priority ($p1).")
	end
	
	# should not allow redispatch from priority: p1 -> p1, or p1 -> p2 -> p1, or p1 -> p2 -> p3 -> p1, etc, otherwise may redispatch indefinitely
	@assert(all(x -> x == 0, redispatch.conditions ^ numPriorities), "Redispatch file has conditions with a cycle that may cause unlimited redispatching.")
	
	return redispatch
end

function readRNetTravelsFile(filename::String)
	rNetTravels = deserializeFile(filename)
	for (i, rNetTravel) in enumerate(rNetTravels)
		@assert(isa(rNetTravel, NetTravel))
		@assert(rNetTravel.isReduced)
		@assert(rNetTravel.modeIndex == i)
	end
	return rNetTravels
end

function readStationsFile(filename::String)
	tables = readTablesFromFile(filename)
	
	table = tables["stations"]
	n = size(table.data,1) # number of stations
	@assert(n >= 1)
	
	# create stations from data in table
	columns = table.columns # shorthand
	stations = Vector{Station}(undef, n)
	for i = 1:n
		stations[i] = Station()
		stations[i].index = columns["index"][i]
		stations[i].location.x = columns["x"][i]
		stations[i].location.y = columns["y"][i]
		stations[i].capacity = columns["capacity"][i]
		
		@assert(stations[i].index == i)
		@assert(stations[i].capacity >= 0)
	end
	
	return stations
end

function readStatsControlFile(filename::String)
	tables = readTablesFromFile(filename)
	
	table = tables["params"]
	param(s::String) = table.columns[s][1]
	warmUpDuration = param("warmUpDuration")
	periodDurationsString = param("periodDurations")
	doCyclePeriodDurations = Bool(param("doCyclePeriodDurations"))
	
	@assert(warmUpDuration >= 0.0)
	
	# create period durations iterator
	@assert(isa(periodDurationsString, String) || isa(periodDurationsString, SubString)) # should be string representation of a vector
	periodDurations = periodDurationsString |> Meta.parse |> eval
	@assert(all(x -> x > 0, periodDurations))
	if doCyclePeriodDurations periodDurations = Iterators.cycle(periodDurations) end
	periodDurationsIter = Stateful(periodDurations)
	
	return periodDurationsIter, warmUpDuration
end

function readTravelFile(filename::String)
	tables = readTablesFromFile(filename)
	
	# create travel object, will populate gradually
	travel = Travel()
	
	# travel modes
	table = tables["travelModes"]
	n = numModes = travel.numModes = size(table.data,1)
	@assert(n >= 1)
	# create travel modes from data in table (network based data will need to be filled later)
	columns = table.columns # shorthand
	travel.modes = Vector{TravelMode}(undef, n)
	for i = 1:n
		travel.modes[i] = TravelMode()
		travel.modes[i].index = columns["travelModeIndex"][i]
		travel.modes[i].offRoadSpeed = columns["offRoadSpeed"][i]
		
		@assert(travel.modes[i].index == i)
		@assert(travel.modes[i].offRoadSpeed > 0)
	end
	
	# travel sets
	table = tables["travelSets"]
	n = size(table.data,1)
	@assert(n >= 1)
	# populate travel object from data in table
	columns = table.columns # shorthand
	numSets = travel.numSets = maximum(columns["travelSetIndex"])
	@assert(sort(unique(columns["travelSetIndex"])) == [1:numSets;]) # all sets from 1:numSets should be used
	@assert(sort(unique(columns["travelModeIndex"])) == [1:numModes;]) # all travel modes should be used
	travel.modeLookup = fill(nullIndex, numSets, numPriorities)
	for i = 1:n
		setIndex = columns["travelSetIndex"][i]
		priority = eval(Meta.parse(columns["priority"][i]))
		modeIndex = columns["travelModeIndex"][i]
		
		@assert(travel.modeLookup[setIndex, Int(priority)] == nullIndex) # value should not yet be filled
		@assert(1 <= modeIndex <= numModes)
		travel.modeLookup[setIndex, Int(priority)] = modeIndex
	end
	@assert(all(travel.modeLookup .!= nullIndex)) # check that travel.modeLookup is filled correctly
	
	# travel sets timing
	table = tables["travelSetsTiming"]
	# populate travel object from data in table
	columns = table.columns # shorthand
	travel.setsStartTimes = columns["startTime"]
	@assert(travel.setsStartTimes[1] > nullTime)
	@assert(issorted(travel.setsStartTimes, lt=<=)) # times should be strictly increasing
	travel.setsTimeOrder = columns["travelSetIndex"]
	@assert(sort(unique(travel.setsTimeOrder)) == [1:numSets;]) # each travel set should be used at least once
	
	# misc
	travel.recentSetsStartTimesIndex = 1
	
	return travel
end

# read deployments from file
function readDeploymentsFile(filename::String)
	tables = readTablesFromFile(filename)
	
	# get counts
	table = tables["miscData"]
	numStations = table.columns["numStations"][1]
	numDeployments = table.columns["numDeployments"][1]
	
	# deployments
	table = tables["deployments"]
	columns = table.columns # shorthand
	# check number of ambulances
	ambIndexCol = columns["ambIndex"]
	@assert(ambIndexCol == collect(1:length(ambIndexCol)))
	# check that table header has all deployments
	for i = 1:numDeployments
		@assert(in("deployment_$i stationIndex", table.header), "Missing deployment $i")
	end
	# create deployments from table
	deployments = Vector{Deployment}(undef, numDeployments)
	for i = 1:numDeployments
		deployments[i] = columns["deployment_$i stationIndex"]
		
		for value in deployments[i]
			@assert(in(value, 1:numStations))
		end
	end
	
	return deployments, numStations
end

function readZhangIpParamsFile(filename::String)
	tables = readTablesFromFile(filename)
	
	zid = ZhangIpData()
	
	table = tables["miscParams"]
	@assert(size(table.data,1) == 1) # table should have one data row
	for fname in [:travelTimeCost, :onRoadMoveUpDiscountFactor, :regretTravelTimeThreshold]
		# setfield!(zid, fname, table.columns[string(fname)][1])
		setfield!(zid, fname, convert(fieldtype(typeof(zid), fname), table.columns[string(fname)][1]))
	end
	if haskey(table.columns, "expectedHospitalHandoverDuration")
		zid.expectedHospitalHandoverDuration = table.columns["expectedHospitalHandoverDuration"][1]
	else # compat
		@assert(haskey(table.columns, "expectedHospitalTransferDuration"))
		@warn("`expectedHospitalTransferDuration` header in Zhang IP params file is deprecated, use `expectedHospitalHandoverDuration` instead.")
		zid.expectedHospitalHandoverDuration = table.columns["expectedHospitalTransferDuration"][1]
	end
	
	table = tables["stationCapacities"]
	n = numStations = size(table.data, 1) # number of stations
	for i = 1:n
		@assert(table.columns["stationIndex"][i] == i)
		push!(zid.stationCapacities, table.columns["capacity"][i])
	end
	
	table = tables["stationMarginalBenefits"]
	n = size(table.data, 1) # number of stations
	@assert(n == numStations)
	for i = 1:n
		@assert(table.columns["stationIndex"][i] == i)
		push!(zid.marginalBenefits, [])
		for j = 1:zid.stationCapacities[i]
			push!(zid.marginalBenefits[i], table.columns["slot_$j"][i])
		end
	end
	
	return zid
end
