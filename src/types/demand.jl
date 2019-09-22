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

# Demand (call) model

# get demand mode for given time and priority
function getDemandMode!(demand::Demand, priority::Priority, startTime::Float)
	@assert(startTime != nullTime)
	@assert(priority != nullPriority)
	
	updateDemandToTime!(demand, startTime) # update demand.recentSetsStartTimesIndex for given start time
	demandSetIndex = demand.setsTimeOrder[demand.recentSetsStartTimesIndex]
	demandModeIndex = demand.modeLookup[demandSetIndex, Int(priority)]
	
	return demand.modes[demandModeIndex]
end

# update demand.recentSetsStartTimesIndex up to given demand start time
function updateDemandToTime!(demand::Demand, startTime::Float)
	# shorthand:
	setsStartTimes = demand.setsStartTimes
	i = demand.recentSetsStartTimesIndex
	n = length(setsStartTimes)
	
	@assert(setsStartTimes[i] <= startTime) # otherwise, have gone back in time?
	
	if i == n || startTime < setsStartTimes[i+1]
		return # do nothing
	end
	
	# use linear search to update position in setsStartTimes
	while i < n && setsStartTimes[i+1] <= startTime
		i += 1
	end
	
	demand.recentSetsStartTimesIndex = i
end

# Creates one set of points to represent the demand in demand.rasters.
# Each point will store one or more demand values, one from each raster.
# Where there are multiple points per raster cell, the demand is divided between the points.
# Points are not created in places where no rasters have demand.
# Requires all rasters to have the same dimensions (only z values can differ).
function createDemandPointsFromRasters(demand::Demand; numCellRows::Int = 1, numCellCols::Int = 1)
	@assert(demand.numRasters >= 1)
	@assert(numCellRows >= 1)
	@assert(numCellCols >= 1)
	
	# shorthand
	@unpack rasters, numRasters = demand
	@unpack nx, ny = rasters[1]
	
	# check that all rasters have same dimensions,
	# otherwise creating points that work for all rasters may be tricky
	for i = 2:numRasters
		@assert(rasters[i].nx == nx)
		@assert(rasters[i].ny == ny)
		@assert(all(j -> isapprox(rasters[i].x[j], rasters[1].x[j]; rtol = eps()), 1:nx))
		@assert(all(j -> isapprox(rasters[i].y[j], rasters[1].y[j]; rtol = eps()), 1:ny))
	end
	
	# create point(s) for each raster cell,
	# except where point is in cell where all rasters have z = 0
	# numCellRows is for y dimension, numCellCols for x
	numCellPoints = numCellRows * numCellCols # number of points per cell
	points = Vector{Point}()
	maxNumPoints = (nx * numCellCols) * (ny * numCellRows)
	rastersPointDemands = [Vector{Float}(undef, maxNumPoints) for i = 1:numRasters] # rastersPointDemands[i][j] gives demand at point j for raster i
	for i = 1:nx, j = 1:ny # raster cell indices
		if all(raster -> raster.z[i,j] == 0, rasters)
			continue
		end
		for cellRow = 1:numCellRows, cellCol = 1:numCellCols
			point = Point()
			point.index = length(points) + 1
			
			# set point location
			point.location.x = rasters[1].x[i] + rasters[1].dx * ((cellCol - 0.5)/numCellCols - 0.5)
			point.location.y = rasters[1].y[j] + rasters[1].dy * ((cellRow - 0.5)/numCellRows - 0.5)
			
			# set demand of point for each raster,
			# divide raster demand by number of points per cell
			for rasterIndex = 1:numRasters
				rastersPointDemands[rasterIndex][point.index] = rasters[rasterIndex].z[i,j] / numCellPoints
			end
			
			push!(points, point)
		end
	end
	
	# remove unneeded values from rastersPointDemands (rastersPointDemands[i][j] where j > length(points))
	numPoints = length(points) # shorthand
	for i = 1:numRasters
		rastersPointDemands[i] = rastersPointDemands[i][1:numPoints]
	end
	
	return points, rastersPointDemands
end

# Sets the nearest node data for the points.
# mutates: points (fields: nearestNodeIndex, nearestNodeDist)
function setPointsNearestNodes!(sim::Simulation, points::Vector{Point})
	# find and set nearest node index of each point
	for point in points
		@assert(point.location.x != nullX && point.location.y != nullY)
		(point.nearestNodeIndex, point.nearestNodeDist) = findNearestNode(sim.map, sim.grid, sim.net.fGraph.nodes, point.location)
	end
end

# Creates an instance of PointsCoverageMode for the given travelMode and coverTime.
function createPointsCoverageMode(sim::Simulation, travelMode::TravelMode, coverTime::Float)
	@assert(coverTime >= 0)
	@assert(!isempty(sim.demandCoverage.nodesPoints))
	
	pointsCoverageMode = PointsCoverageMode()
	pointsCoverageMode.coverTime = coverTime
	pointsCoverageMode.travelMode = travelMode
	pointsCoverageMode.points = sim.demandCoverage.points
	
	# shorthand:
	points = pointsCoverageMode.points
	numPoints = length(points)
	@unpack stations, numStations = sim
	numNodes = length(sim.net.fGraph.nodes)
	nodesPoints = sim.demandCoverage.nodesPoints
	
	# # determine point coverage
	# stationsCoverPoints = fill(false, numStations, numPoints)
	# for i = 1:numStations
		# station = stations[i] # shorthand
		
		# # get travel time from station to nearest node
		# (node1, dist1) = (station.nearestNodeIndex, station.nearestNodeDist)
		# time1 = offRoadTravelTime(travelMode, dist1) # time to reach nearest node
		
		# for j = 1:numPoints
			# # get shortest path travel time
			# point = points[j] # shorthand
			# (node2, dist2) = (point.nearestNodeIndex, point.nearestNodeDist)
			# pathTravelTime = shortestPathTravelTime(net, travelMode.index, node1, node2)
			# time2 = offRoadTravelTime(travelMode, dist2) # time to reach nearest node
			# if time1 + pathTravelTime + time2 <= coverTime
				# stationsCoverPoints[i,j] = true
			# end
		# end
	# end
	
	# determine point coverage
	stationsCoverPoints = fill(false, numStations, numPoints)
	for i = 1:numStations
		station = stations[i] # shorthand
		
		# get travel time from station to nearest node
		(node1, dist1) = (station.nearestNodeIndex, station.nearestNodeDist)
		time1 = offRoadTravelTime(travelMode, dist1) # time to reach nearest node from station
		
		for j = 1:numNodes
			pathTravelTime = shortestPathTravelTime(sim.net, travelMode.index, node1, j)
			if time1 + pathTravelTime <= coverTime
				for k in nodesPoints[j] # indices of points for which node j is the nearest node
					point = points[k] # shorthand
					@assert(point.nearestNodeIndex == j)
					time2 = offRoadTravelTime(travelMode, point.nearestNodeDist) # time to reach point from nearest node
					stationsCoverPoints[i,k] = (time1 + pathTravelTime + time2 <= coverTime)
				end
			end
		end
	end
	
	# group points with coverage by same set of stations
	stationSetsPoints = Dict{Vector{Bool},Vector{Int}}() # stationSetsPoints[stationSet] gives the points that belong to the stationSet
	for j = 1:numPoints
		stationSet = stationsCoverPoints[:,j] # stationSet[i] = if station i covers point j
		# will put points that are not covered into a set, though alternatively they could be ignored
		stationSets = get!(stationSetsPoints, stationSet, Int[])
		push!(stationSets, j)
	end
	pointSets = collect(values(stationSetsPoints)) # pointSets[i] has set of all point indices covered by the same unique set of stations
	stationSets = [findall(stationSet) for stationSet in keys(stationSetsPoints)] # stationSets[i] = unique set of stations for pointSets[i]
	# according to documentation, 'values' and 'keys' return elements in the same order, so pointSets[i] corresponds with stationSets[i]
	
	stationsCoverPointSets = Vector{Vector{Int}}(undef, numStations) # stationsCoverPointSets[i] = list of point set indices covered by station i
	for i = 1:numStations
		stationsCoverPointSets[i] = findall(stationSet -> in(i,stationSet), stationSets)
	end
	
	pointsCoverageMode.pointSets = pointSets
	pointsCoverageMode.stationSets = stationSets
	pointsCoverageMode.stationsCoverPointSets = stationsCoverPointSets
	
	return pointsCoverageMode # will still need to set pointsCoverageMode.index
end

function getPointsCoverageMode(demandCoverage::DemandCoverage, travelMode::TravelMode, coverTime::Float)
	dc = demandCoverage # shorthand
	@assert(travelMode.index <= length(dc.pointsCoverageModeLookup))
	@assert(haskey(dc.pointsCoverageModeLookup[travelMode.index], coverTime))
	pointsCoverageMode = dc.pointsCoverageModes[dc.pointsCoverageModeLookup[travelMode.index][coverTime]]
	return pointsCoverageMode
end

# mutates: sim.travel state
function getPointsCoverageMode!(sim::Simulation, demandPriority::Priority, currentTime::Float)
	@assert(demandPriority != nullPriority)
	@assert(currentTime != nullTime)
	travelPriority = sim.responseTravelPriorities[demandPriority]
	travelMode = getTravelMode!(sim.travel, travelPriority, currentTime)
	coverTime = sim.demandCoverage.coverTimes[demandPriority]
	pointsCoverageMode = getPointsCoverageMode(sim.demandCoverage, travelMode, coverTime)
	@assert(pointsCoverageMode != nothing)
	return pointsCoverageMode
end

# mutates: sim.travel and sim.demand states
function getPointSetsDemands!(sim::Simulation, demandPriority::Priority, currentTime::Float;
	pointsCoverageMode::Union{PointsCoverageMode,Nothing} = nothing)
	if pointsCoverageMode == nothing
		pointsCoverageMode = getPointsCoverageMode!(sim, demandPriority, currentTime)
	end
	demandMode = getDemandMode!(sim.demand, demandPriority, currentTime)
	pointSetsDemands = sim.demandCoverage.pointSetsDemands[pointsCoverageMode.index, demandMode.rasterIndex]
	return pointSetsDemands
end

"""
	function initDemand!(sim::Simulation, demand::Union{Demand,Nothing} = nothing;
		demandFilename::String = "")
Initialise demand data for `sim`.
Will initialise from the first item in this list: `demand`, `demandFilename` file, `sim.demand` (if set), `sim.inputFiles[\"demand\"].path` file.
Deletes old demand coverage data in `sim.demandCoverage`.

Mutates: `sim.demand`, `sim.demandCoverage`
"""
function initDemand!(sim::Simulation, demand::Union{Demand,Nothing} = nothing;
	demandFilename::String = "")
	if demand == nothing
		if demandFilename != ""
			demand = readDemandFile(demandFilename)
		elseif sim.demand.numSets == nullIndex
			if haskey(sim.inputFiles, "demand")
				demand = readDemandFile(sim.inputFiles["demand"].path)
			else
				error("no demand data given")
			end
		else
			demand = sim.demand
		end
	end
	demand.initialised = true
	sim.demand = demand
	
	# old demand coverage data may no longer be valid
	sim.demandCoverage = DemandCoverage(sim.demandCoverage) # keep some old params
end

# Initialises the points coverage modes data in sim.demandCoverage.
# mutates: sim.demandCoverage (fields: pointsCoverageModes, pointsCoverageModeLookup)
function initPointsCoverageModes!(sim::Simulation)
	# shorthand:
	dc = sim.demandCoverage
	travel = sim.travel
	
	# go through all travel sets and demand priorities, create points coverage mode for each call response travel priority and cover time
	# dc.pointsCoverageModes
	dc.pointsCoverageModeLookup = [Dict{Float,Int}() for i = 1:travel.numModes]
	for i = 1:travel.numSets, demandPriority in priorities
		# get travel mode and cover time
		travelPriority = sim.responseTravelPriorities[demandPriority]
		travelMode = travel.modes[travel.modeLookup[i,Int(travelPriority)]]
		coverTime = sim.demandCoverage.coverTimes[demandPriority]
		
		# if combination of travel mode and cover time is new, need to calculate points coverage
		if !haskey(dc.pointsCoverageModeLookup[travelMode.index], coverTime)
			pointsCoverageMode = createPointsCoverageMode(sim, travelMode, coverTime)
			pointsCoverageMode.index = length(dc.pointsCoverageModes) + 1
			push!(dc.pointsCoverageModes, pointsCoverageMode)
			dc.pointsCoverageModeLookup[travelMode.index][coverTime] = pointsCoverageMode.index
		end
	end
end

"""
	function initDemandCoverage!(sim::Simulation;
		coverTimes::Union{Dict{Priority,Float},Nothing} = nothing,
		rasterCellNumRows::Int = 1, rasterCellNumCols::Int = 1)
Initialises demand coverage data for `sim`.
If `coverTimes` is given, will use this (along with `rasterCellNumRows` and `rasterCellNumCols`), otherwise will use the previous parameter values set in `sim.demandCoverage`.

# Keyword arguments
- `coverTimes` is the target coverage time for each demand priority, e.g. `coverTimes = Dict([p => 8/60/24 for p in priorities])` sets a target time of 8 minutes (converted to days) for each priority. Can omit this if `sim.demandCoverage.coverTimes` is already set.
- `rasterCellNumRows` and `rasterCellNumCols` are used to set the number of demand points to represent demand per raster cell for `sim.demand.rasters`

Mutates: `sim.demandCoverage`
"""
function initDemandCoverage!(sim::Simulation;
	coverTimes::Union{Dict{Priority,Float},Nothing} = nothing,
	rasterCellNumRows::Int = 1, rasterCellNumCols::Int = 1)
	
	@assert(!sim.used)
	@assert(sim.travel.recentSetsStartTimesIndex == 1)
	sim.demand.initialised || initDemand!(sim)
	@assert(sim.demand.recentSetsStartTimesIndex == 1)
	
	# shorthand
	@unpack demand, travel = sim
	@unpack rasters, numRasters = demand
	numNodes = length(sim.net.fGraph.nodes)
	
	# set demand coverage params
	oldDemandCoverage = sim.demandCoverage
	if coverTimes != nothing
		sim.demandCoverage = DemandCoverage(coverTimes, rasterCellNumRows, rasterCellNumCols)
	else # use old params
		@assert(!isempty(oldDemandCoverage.coverTimes))
		sim.demandCoverage = DemandCoverage(oldDemandCoverage)
	end
	dc = sim.demandCoverage # shorthand
	demandPriorities = priorities
	@assert(all(p -> haskey(dc.coverTimes, p), demandPriorities), "coverTimes not set")
	@assert(dc.rasterCellNumRows >= 1)
	@assert(dc.rasterCellNumCols >= 1)
	
	# create demand points
	dc.points, dc.rastersPointDemands = createDemandPointsFromRasters(sim.demand; numCellRows = dc.rasterCellNumRows, numCellCols = dc.rasterCellNumCols)
	points = dc.points # shorthand
	setPointsNearestNodes!(sim, points)
	numPoints = length(points) # shorthand
	
	# set nodesPoints
	dc.nodesPoints = [Int[] for i = 1:numNodes]
	for point in points
		push!(dc.nodesPoints[point.nearestNodeIndex], point.index)
	end
	
	# set points coverage modes
	initPointsCoverageModes!(sim)
	numPointsCoverageModes = length(dc.pointsCoverageModes)
	
	# set pointSetsDemands
	# note that pointSetsDemands does not need to be set for all combinations of pointsCoverageMode.index and rasterIndex, as some combinations may not occur in simulation
	dc.pointSetsDemands = [Vector{Float}() for i = 1:numPointsCoverageModes, j = 1:demand.numRasters]
	times = vcat(travel.setsStartTimes, demand.setsStartTimes) |> unique |> sort
	for time in times, demandPriority in demandPriorities
		pointsCoverageMode = getPointsCoverageMode!(sim, demandPriority, time)
		demandMode = getDemandMode!(demand, demandPriority, time)
		rasterIndex = demandMode.rasterIndex # shorthand
		
		# if combination of pointsCoverageMode and rasterIndex is new, need to calculate point sets demands
		pointDemands = dc.rastersPointDemands[rasterIndex]
		if isempty(dc.pointSetsDemands[pointsCoverageMode.index, rasterIndex])
			pointSetsDemand = zeros(Float, length(pointsCoverageMode.pointSets))
			for (pointSetIndex, pointSet) in enumerate(pointsCoverageMode.pointSets)
				for pointIndex in pointSet
					pointSetsDemand[pointSetIndex] += pointDemands[pointIndex]
				end
			end
			dc.pointSetsDemands[pointsCoverageMode.index, rasterIndex] = pointSetsDemand
		end
	end
	# reset travel and demand state (previously changed by getTravelMode! and getDemandMode!)
	travel.recentSetsStartTimesIndex = 1
	demand.recentSetsStartTimesIndex = 1
	
	dc.initialised = true
end

# Calculate, for a given points coverage mode, and number of ambulances assigned to each
# station, how many ambulances cover each point set.
function calcPointSetsCoverCounts(pointsCoverageMode::PointsCoverageMode, stationsNumAmbs::Vector{Int})
	# stationsNumAmbs[i] gives the number of ambulances that can provide coverage and are assigned to station i
	
	@assert(all(n -> n >= 0, stationsNumAmbs))
	
	# count how many times each point set is covered
	numPointSets = length(pointsCoverageMode.pointSets) # shorthand
	pointSetsCoverCount = zeros(Int, numPointSets)
	for i = 1:numPointSets
		for j in pointsCoverageMode.stationSets[i]
			pointSetsCoverCount[i] += stationsNumAmbs[j]
		end
	end
	
	# # alternative code:
	# @assert(length(stationsNumAmbs) == length(pointsCoverageMode.stationsCoverPointSets))
	# pointSetsCoverCount = zeros(Int, length(pointsCoverageMode.pointSets))
	# for (i,stationNumAmbs) in enumerate(stationsNumAmbs)
		# for j in pointsCoverageMode.stationsCoverPointSets[i]
			# pointSetsCoverCount[j] += stationNumAmbs
		# end
	# end
	
	return pointSetsCoverCount
end

# Calculate, for each demand priority, how many ambulances cover each point set in
# pointsCoverageMode.pointSet, where pointsCoverageMode is for the demand priority and current time.
# mutates: sim.travel state
function calcPointSetsCoverCounts!(sim::Simulation, currentTime::Float, stationsNumAmbs::Vector{Int};
	demandPriorities::Vector{Priority} = priorities)
	# stationsNumAmbs[i] gives the number of ambulances that can provide coverage and are assigned to station i
	
	@assert(length(stationsNumAmbs) == sim.numStations)
	@assert(all(n -> n >= 0, stationsNumAmbs))
	
	demandsPointSetsCoverCounts = Dict{Priority,Vector{Int}}() # demandsPointSetsCoverCounts[p][i] gives the number of ambulances assigned to stations that cover point set i for demand priority p
	pointsCoverageModesUsed = Dict{PointsCoverageMode,Priority}() # keep track of which points coverage mode are used, and any demand priority that uses each
	for demandPriority in demandPriorities
		pointsCoverageMode = getPointsCoverageMode!(sim, demandPriority, currentTime)
		if haskey(pointsCoverageModesUsed, pointsCoverageMode)
			# have already calculated the point sets cover counts for this points coverage mode
			priority = pointsCoverageModesUsed[pointsCoverageMode]
			demandsPointSetsCoverCounts[demandPriority] = demandsPointSetsCoverCounts[priority] # should use copy(), so that vectors are independent?
		else
			demandsPointSetsCoverCounts[demandPriority] = calcPointSetsCoverCounts(pointsCoverageMode, stationsNumAmbs)
			pointsCoverageModesUsed[pointsCoverageMode] = demandPriority
		end
	end
	
	return demandsPointSetsCoverCounts
end

# Calculates, for each demand priority, how much of the demand is covered by 1,2,... ambulances
# mutates: sim.travel and sim.demand states
function calcDemandCoverCounts!(sim::Simulation, currentTime::Float, stationsNumAmbs::Vector{Int};
	demandPriorities::Vector{Priority} = priorities)
	# stationsNumAmbs[i] gives the number of ambulances that can provide coverage and are assigned to station i
	
	@assert(length(stationsNumAmbs) == sim.numStations)
	@assert(all(n -> n >= 0, stationsNumAmbs))
	
	# shorthand
	numAmbs = sum(stationsNumAmbs)
	
	demandsPointSetsCoverCounts = calcPointSetsCoverCounts!(sim, currentTime, stationsNumAmbs; demandPriorities = demandPriorities)
	
	demandCoverCounts = Dict{Priority,Vector{Float}}() # demandCoverCounts[p][i] gives the amount of demand covered by i ambulances for demand priority p
	coverCountsCalculated = [Dict{Vector{Int},Vector{Int}}() for i = 1:sim.demand.numRasters] # coverCountsCalculated[rasterIndex][pointSetsCoverCount] stores the coverCounts (before multiplying by demandMode.rasterMultiplier) for the combination of args
	for demandPriority in demandPriorities
		pointsCoverageMode = getPointsCoverageMode!(sim, demandPriority, currentTime)
		demandMode = getDemandMode!(sim.demand, demandPriority, currentTime)
		rasterIndex = demandMode.rasterIndex # shorthand
		
		pointSetsCoverCount = demandsPointSetsCoverCounts[demandPriority]
		pointSetsDemands = sim.demandCoverage.pointSetsDemands[pointsCoverageMode.index, rasterIndex]
		
		coverCounts = nothing # init
		if haskey(coverCountsCalculated[rasterIndex], pointSetsCoverCount)
			coverCounts = coverCountsCalculated[rasterIndex][pointSetsCoverCount]
		else
			coverCounts = zeros(Float, numAmbs) # coverCounts[i] is the demand covered by i ambulances
			for (i,pointSet) in enumerate(pointsCoverageMode.pointSets)
				if pointSetsCoverCount[i] > 0
					coverCounts[pointSetsCoverCount[i]] += pointSetsDemands[i]
				end
			end
			coverCountsCalculated[rasterIndex][pointSetsCoverCount] = coverCounts
		end
		demandCoverCounts[demandPriority] = coverCounts * demandMode.rasterMultiplier
	end
	
	return demandCoverCounts
end
