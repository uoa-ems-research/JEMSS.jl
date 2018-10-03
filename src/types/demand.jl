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
	rasters = demand.rasters
	numRasters = demand.numRasters
	nx = rasters[1].nx
	ny = rasters[1].ny
	
	# check that all rasters have same dimensions,
	# otherwise creating points that work for all rasters may be tricky
	for i = 2:numRasters
		@assert(rasters[i].nx == nx)
		@assert(rasters[i].ny == ny)
		@assert(all(j -> isapprox(rasters[i].x[j], rasters[1].x[j]; rtol = eps()), j = 1:nx))
		@assert(all(j -> isapprox(rasters[i].y[j], rasters[1].y[j]; rtol = eps()), j = 1:ny))
	end
	
	# create point(s) for each raster cell,
	# except where point is in cell where all rasters have z = 0
	# numCellRows is for y dimension, numCellCols for x
	numCellPoints = numCellRows * numCellCols # number of points per cell
	points = Vector{Point}()
	for i = 1:nx, j = 1:ny # raster cell indices
		if all(raster -> raster.z[i,j] == 0, rasters)
			continue
		end
		for cellRow = 1:numCellRows, cellCol = 1:numCellCols
			point = Point()
			point.index = length(points) + 1
			
			# set point location
			x = rasters[1].x[i] + rasters[1].dx * ((cellCol - 0.5)/numCellCols - 0.5)
			y = rasters[1].y[j] + rasters[1].dy * ((cellRow - 0.5)/numCellRows - 0.5)
			point.location = Location(x,y)
			
			# set demand of point for each raster,
			# divide raster demand by number of points per cell
			point.value = Dict{String,Any}()
			point.value["demands"] = Vector{Float}(numRasters)
			for rasterIndex = 1:numRasters
				point.value["demands"][rasterIndex] = rasters[rasterIndex].z[i,j] / numCellPoints
			end
			
			push!(points, point)
		end
	end
	
	return points
end

# Sets the nearest node data for the points.
# mutates: points (fields: nearestNodeIndex, nearestNodeDist)
function setPointsNearestNodes!(sim::Simulation, points::Vector{Point})
	# find and set nearest node index of each point
	for point in points
		@assert(point.location.x != nullX && point.location.y != nullY)
		(point.nearestNodeIndex, point.nearestNodeDist) = findNearestNodeInGrid(sim.map, sim.grid, sim.net.fGraph.nodes, point.location)
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
	stations = sim.stations
	numStations = length(stations)
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
			# (pathTravelTime, rNodes) = shortestPathTravelTime(net, travelMode.index, node1, node2)
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
			(pathTravelTime, rNodes) = shortestPathTravelTime(sim.net, travelMode.index, node1, j)
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
	
	# to do: speed up grouping of points into sets (below)
	
	# group points with coverage by same set of stations
	pointSets = Vector{Set{Int}}() # pointSets[i] has set of all point indices covered by the same unique set of stations
	stationSets = Vector{Set{Int}}() # stationSets[i] = unique set of stations for pointSets[i]
	for j = 1:numPoints
		stationsCoverPoint = Set(find(stationsCoverPoints[:,j])) # indices of stations that cover point j
		# if !isempty(stationsCoverPoint) # point is not covered, ignore it?
		if !in(stationsCoverPoint, stationSets)
			push!(stationSets, stationsCoverPoint)
			push!(pointSets, Set())
		end
		# add point to pointSets for corresponding stationSets
		k = findfirst(stationSets, stationsCoverPoint)
		push!(pointSets[k], j)
		# end
	end
	
	stationsCoverPointSets = Vector{Vector{Int}}(numStations) # stationsCoverPointSets[i] = list of point set indices covered by station i
	for i = 1:numStations
		stationsCoverPointSets[i] = find(stationSet -> in(i,stationSet), stationSets)
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
	coverTime = sim.demandCoverTimes[demandPriority]
	pointsCoverageMode = getPointsCoverageMode(sim.demandCoverage, travelMode, coverTime)
	@assert(pointsCoverageMode != nothing)
	return pointsCoverageMode
end

# mutates: sim.travel and sim.demand states
function getPointSetsDemands!(sim::Simulation, demandPriority::Priority, currentTime::Float;
	pointsCoverageMode::Union{PointsCoverageMode,Void} = nothing)
	if pointsCoverageMode == nothing
		pointsCoverageMode = getPointsCoverageMode!(sim, demandPriority, currentTime)
	end
	demandMode = getDemandMode!(sim.demand, demandPriority, currentTime)
	pointSetsDemands = sim.demandCoverage.pointSetsDemands[pointsCoverageMode.index, demandMode.rasterIndex]
	return pointSetsDemands
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
	priorities = setdiff([instances(Priority)...], [nullPriority])
	for i = 1:travel.numSets, demandPriority in priorities
		# get travel mode and cover time
		travelPriority = sim.responseTravelPriorities[demandPriority]
		travelMode = travel.modes[travel.modeLookup[i,Int(travelPriority)]]
		coverTime = sim.demandCoverTimes[demandPriority]
		
		# if combination of travel mode and cover time is new, need to calculate points coverage
		if !haskey(dc.pointsCoverageModeLookup[travelMode.index], coverTime)
			pointsCoverageMode = createPointsCoverageMode(sim, travelMode, coverTime)
			pointsCoverageMode.index = length(dc.pointsCoverageModes) + 1
			push!(dc.pointsCoverageModes, pointsCoverageMode)
			dc.pointsCoverageModeLookup[travelMode.index][coverTime] = pointsCoverageMode.index
		end
	end
end

# Initialises data in sim.demandCoverage.
# mutates: sim.demandCoverage
function initDemandCoverage!(sim::Simulation;
	rasterCellNumRows::Int = 1, rasterCellNumCols::Int = 1)
	@assert(!sim.used)
	@assert(sim.travel.recentSetsStartTimesIndex == 1)
	@assert(sim.demand.recentSetsStartTimesIndex == 1)
	
	# shorthand
	demand = sim.demand
	rasters = demand.rasters
	numRasters = demand.numRasters
	numNodes = length(sim.net.fGraph.nodes)
	travel = sim.travel
	
	# create demand points
	dc = demandCoverage = sim.demandCoverage = DemandCoverage()
	points = dc.points = createDemandPointsFromRasters(sim.demand; numCellRows = rasterCellNumRows, numCellCols = rasterCellNumCols) # should add kwargs
	setPointsNearestNodes!(sim, points)
	numPoints = length(points) # shorthand
	
	# set nodesPoints
	dc.nodesPoints = [Int[] for i = 1:numNodes]
	for point in points
		push!(dc.nodesPoints[point.nearestNodeIndex], point.index)
	end
	
	# set rastersPointDemands
	dc.rastersPointDemands = [Vector{Float}(numPoints) for i = 1:numRasters]
	for i = 1:numRasters, j = 1:numPoints
		dc.rastersPointDemands[i][j] = points[j].value["demands"][i]
	end
	
	# set points coverage modes
	initPointsCoverageModes!(sim)
	numPointsCoverageModes = length(dc.pointsCoverageModes)
	
	# set pointSetsDemands
	# note that pointSetsDemands does not need to be set for all combinations of pointsCoverageMode.index and rasterIndex, as some combinations may not occur in simulation
	dc.pointSetsDemands = [Vector{Float}() for i = 1:numPointsCoverageModes, j = 1:demand.numRasters]
	times = vcat(travel.setsStartTimes, demand.setsStartTimes) |> unique |> sort
	priorities = setdiff([instances(Priority)...], [nullPriority])
	for time in times, demandPriority in priorities
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
end

# Calculate, for a given points coverage mode, and number of available ambulances assigned to each
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
	demandPriorities::Vector{Priority} = setdiff([instances(Priority)...], [nullPriority]))
	# stationsNumAmbs[i] gives the number of ambulances that can provide coverage and are assigned to station i
	
	@assert(length(stationsNumAmbs) == length(sim.stations))
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
	demandPriorities::Vector{Priority} = setdiff([instances(Priority)...], [nullPriority]))
	# stationsNumAmbs[i] gives the number of ambulances that can provide coverage and are assigned to station i
	
	@assert(length(stationsNumAmbs) == length(sim.stations))
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
