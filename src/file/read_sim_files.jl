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
	ambulances = Vector{Ambulance}(n)
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
	fieldNames = setdiff(table.header, keepAllFields ? [] : ["index", "fromNode", "toNode", [string("mode_", i) for i = 1:numModes]...])
	rowsFields = tableRowsFieldDicts(table, fieldNames)
	# get arc data from table
	for i = 1:numArcs
		arcs[i] = Arc()
		arcs[i].index = indexCol[i]
		arcs[i].fromNodeIndex = fromNodeCol[i]
		arcs[i].toNodeIndex = toNodeCol[i]
		arcs[i].fields = rowsFields[i]
		# read travel times
		for j = 1:numModes
			travelTimes[j,i] = modeCols[j][i]
		end
		@assert(arcs[i].index == i)
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
			arcs[j].fields = deepcopy(arcs[i].fields)
		end
		# copy travel times
		travelTimes[:, (1:numArcs) + numArcs] = travelTimes[:, 1:numArcs]
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
		
		@assert(calls[i].index == i)
		@assert(calls[i].priority != nullPriority)
		@assert(calls[i].arrivalTime >= 0)
		@assert(calls[i].dispatchDelay >= 0)
		@assert(calls[i].onSceneDuration >= 0)
		if calls[i].transfer
			@assert(calls[i].transferDuration >= 0)
		end
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
	compTable = CompTable(numAmbs, numStations)
	for i in columns["numAmbs"]
		for j = 1:numStations
			compTable[i,j] = columns[string("station_", j)][i]
		end
	end
	
	@assert(all(compTable .>= 0))
	
	# check that compliance table rows add up to row indices
	@assert(vec(sum(compTable, 2)) == [1:numAmbs;])
	
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
	demand.rasters = Vector{Raster}(n)
	demand.rasterFilenames = Vector{String}(n)
	@assert(allunique(columns["rasterFilename"]), "There are duplicate raster filenames in the demand file, this is unnecessary.")
	for i = 1:n
		@assert(columns["rasterIndex"][i] == i)
		
		rasterFilename = String(columns["rasterFilename"][i])
		@assert(isfile(rasterFilename))
		demand.rasterFilenames[i] = rasterFilename
		demand.rasters[i] = readRasterFile(rasterFilename)
	end
	
	# demand modes
	table = tables["demandModes"]
	n = numModes = demand.numModes = size(table.data,1)
	@assert(n >= 1)
	columns = table.columns # shorthand
	demand.modes = Vector{DemandMode}(n)
	for i = 1:n
		@assert(columns["modeIndex"][i] == i)
		
		demand.modes[i] = DemandMode()
		demand.modes[i].index = i
		rasterIndex = demand.modes[i].rasterIndex = columns["rasterIndex"][i]
		@assert(1 <= rasterIndex <= numRasters)
		demand.modes[i].raster = demand.rasters[rasterIndex]
		demand.modes[i].priority = columns["priority"][i] |> parse |> eval
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
		modeIndices = columns["modeIndices"][i] |> parse |> eval
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
		priority = eval(parse(columns["demandPriority"][i]))
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
	
	# event dictionary
	table = tables["eventDict"]
	eventNames = table.columns["name"]
	eventKeys = table.columns["key"]
	n = length(eventKeys) # number of event types
	@assert(length(eventNames) == n)
	@assert(allunique(eventNames))
	@assert(allunique(eventKeys))
	eventDict = Dict{Int,EventForm}()
	for i = 1:n
		eventDict[eventKeys[i]] = eval(parse(eventNames[i]))
	end
	
	# create events from data in events table
	table = tables["events"]
	c = table.columns # shorthand
	(eventIndexCol = c["index"]); (parentIndexCol = c["parentIndex"]); (eventKeyCol = c["eventKey"]); (timeCol = c["time"]); (ambIndexCol = c["ambIndex"]); (callIndexCol = c["callIndex"]); (stationIndexCol = c["stationIndex"]) # shorthand, to avoid repeated dict lookups
	data = table.data # shorthand
	fileEnded = (data[end,1] == "end") # true if writing all events to file ended before file was closed
	n = size(data,1) - fileEnded # number of events
	events = Vector{Event}(n)
	for i = 1:n
		events[i] = Event()
		events[i].index = eventIndexCol[i]
		events[i].parentIndex = parentIndexCol[i]
		events[i].form = eventDict[eventKeyCol[i]]
		events[i].time = timeCol[i]
		events[i].ambIndex = ambIndexCol[i]
		events[i].callIndex = callIndexCol[i]
		events[i].stationIndex = stationIndexCol[i]
		
		@assert(events[i].form != nullEvent)
		@assert(events[i].time > nullTime)
		@assert(events[i].ambIndex != nullIndex || events[i].callIndex != nullIndex)
	end
	
	# set eventsChildren; eventsChildren[i] gives events that are children of the ith event that occurred
	eventsChildren = [Vector{Event}() for i = 1:n]
	for i = 1:n
		j = events[i].parentIndex
		if j != nullIndex
			push!(eventsChildren[j], events[i])
		end
	end
	
	return events, eventsChildren, fileEnded, inputFiles, fileChecksums
end

# read geographic data from file and apply f to data
function readGeoFile(f::Function, filename::String)
	return ArchGDAL.registerdrivers() do
		ArchGDAL.read(filename) do dataset
			f(dataset)
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
	hospitals = Vector{Hospital}(n)
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
	nodes = Vector{Node}(n)
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
	targetResponseTimes = Vector{Float}(n)
	responseTravelPriorities = Dict([p => p for p in priorities]) # default is to have response travel priority equal to call priority
	for i = 1:n
		@assert(columns["priority"][i] == i)
		@assert(eval(parse(columns["name"][i])) == Priority(i))
		targetResponseTimes[i] = columns["targetResponseTime"][i]
		
		if haskey(columns, "responseTravelPriority")
			responseTravelPriorities[Priority(i)] = eval(parse(columns["responseTravelPriority"][i]))
		end
	end
	
	# target response times should be positive
	@assert(all(targetResponseTimes .> 0))
	
	return targetResponseTimes, responseTravelPriorities
end

function readPriorityListFile(filename::String)
	tables = readTablesFromFile(filename)
	
	table = tables["priorityList"]
	n = size(table.data,1) # number of ambulances
	@assert(n >= 1)
	
	# create priority list from data in table
	columns = table.columns # shorthand
	priorityList = Vector{Int}(n)
	for i = 1:n
		@assert(columns["numAmbs"][i] == i)
		priorityList[i] = columns["stationIndex"][i]
	end
	
	# station indices should be positive
	@assert(all(priorityList .>= 1))
	
	return priorityList
end

# read raster file using ArchGDAL package, return as custom Raster type
function readRasterFile(rasterFilename::String)
	
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
	# (xMin, dx, z) = (dx > 0) ? (x1 + 0.5*dx, dx, z) : (x1 + (nx-0.5)*dx, -dx, flipdim(z,1))
	(yMin, dy, z) = (dy > 0) ? (y1 + 0.5*dy, dy, z) : (y1 + (ny-0.5)*dy, -dy, flipdim(z,2))
	x = collect(range(xMin, dx, nx))
	y = collect(range(yMin, dy, ny))
	
	return Raster(x, y, z)
end

function readRNetTravelsFile(filename::String)
	rNetTravels = deserializeFile(filename)
	for rNetTravel in rNetTravels
		@assert(isa(rNetTravel, NetTravel))
		@assert(rNetTravel.isReduced)
		@assert(isa(rNetTravel.spTimes, Array{FloatSpTime,2}))
		@assert(isa(rNetTravel.spFadjIndex, Array{IntFadj,2}))
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
	stations = Vector{Station}(n)
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
	travel.modes = Vector{TravelMode}(n)
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
		priority = eval(parse(columns["priority"][i]))
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
	deployments = Vector{Deployment}(numDeployments)
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
	for fname in [:travelTimeCost, :onRoadMoveUpDiscountFactor, :regretTravelTimeThreshold, :expectedHospitalTransferDuration]
		# setfield!(zid, fname, table.columns[string(fname)][1])
		setfield!(zid, fname, convert(fieldtype(typeof(zid), fname), table.columns[string(fname)][1]))
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
