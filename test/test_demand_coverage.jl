# create test data
info("Creating test data")
testRegionDataFolder = "data/regions/small/1"
runGenConfig(joinpath(testRegionDataFolder, "gen_config.xml"), overwriteOutputPath = true, doPrint = false)

# initialise sim
info("Initialising sim")
sim = initSimulation(joinpath(testRegionDataFolder, "sim_config.xml"), doPrint = false);

# initialise demand and demand coverage
sim.demand = readDemandFile(joinpath(testRegionDataFolder, "demand", "demand.csv"))
sim.demandCoverTimes = Dict([p => 20/60/24 for p in instances(Priority)])
sim.demandCoverTimes[highPriority] = 12/60/24
JEMSS.initDemandCoverage!(sim, rasterCellNumRows = 2, rasterCellNumCols = 2)

# shorthand
demand = sim.demand
demandCoverage = sim.demandCoverage;
nodesPoints = demandCoverage.nodesPoints
points = demandCoverage.points
numPoints = length(points)
nodes = sim.net.fGraph.nodes
numNodes = length(nodes)
pointsCoverageModes = demandCoverage.pointsCoverageModes;
stations = sim.stations
numStations = length(sim.stations)
pointsCoverageModeLookup = demandCoverage.pointsCoverageModeLookup
travel = sim.travel;
priorities = setdiff([instances(Priority)...], [nullPriority])

@testset "demand coverage init" begin
	# check sim.demandCoverage after being initialised with function initDemandCoverage!
	
	# check that demand of points adds up to demand of raster, for each raster
	# note that this is not the true demand, as it has not been multiplied by a DemandMode.rasterMultiplier value
	for i = 1:demand.numRasters
		@test isapprox(sum(demand.rasters[i].z), sum(demandCoverage.rastersPointDemands[i]), rtol = eps(Float))
	end
	
	# check demandCoverage.nodesPoints
	@assert(all(i -> points[i].index == i, 1:numPoints))
	@test vcat(nodesPoints...) |> sort == [1:numPoints;] # each point is included exactly once
	@test all(i -> all(j -> points[j].nearestNodeIndex == i, nodesPoints[i]), 1:numNodes) # nodesPoints matches data stored in point.nearestNodeIndex for all points
	
	# check demandCoverage.pointsCoverageModes
	@test allunique([[pcm.coverTime, pcm.travelMode.index] for pcm in pointsCoverageModes]) # each mode should have a unique combination of coverTime and travelMode
	for (i,pointsCoverageMode) in enumerate(pointsCoverageModes)
		@assert pointsCoverageMode.index == i
		@assert pointsCoverageMode.points == demandCoverage.points
		
		# will check fields of pointsCoverageMode: pointSets, stationSets, stationsCoverPointSets
		
		# shorthand
		pointSets = pointsCoverageMode.pointSets
		stationSets = pointsCoverageMode.stationSets
		stationsCoverPointSets = pointsCoverageMode.stationsCoverPointSets
		
		@assert(length(pointSets) == length(stationSets))
		@test vcat(pointSets...) |> sort == [1:numPoints;] # each point should be included in exactly one point set
		@test allunique(stationSets) # each set of stations should be unique, otherwise any duplicates should be combined
		
		travelMode = pointsCoverageMode.travelMode # shorthand
		stationsCoverPoints = fill(false, numStations, numPoints) # stationsCoverPoints[i,j] is true if station i covers point j
		@assert(all(i -> stations[i].index == i, 1:numStations))
		for (i,station) in enumerate(stations)
			for (j,point) in enumerate(points)
				# calculate travel time from station to point
				time1 = offRoadTravelTime(travelMode, station.nearestNodeDist)
				pathTravelTime = shortestPathTravelTime(sim.net, travelMode.index, station.nearestNodeIndex, point.nearestNodeIndex)
				time2 = offRoadTravelTime(travelMode, point.nearestNodeDist)
				travelTime = time1 + pathTravelTime + time2
				
				# check if station covers point
				if travelTime <= pointsCoverageMode.coverTime
					stationsCoverPoints[i,j] = true
				end
			end
		end
		
		@test all(i -> all(j -> find(stationsCoverPoints[:,j]) == stationSets[i], pointSets[i]), 1:length(pointSets)) # assume that both vectors (stationsCoverPoints[:,j] and stationSets[i]) are sorted the same
	end
	
	# check pointsCoverageModeLookup
	for pcm in pointsCoverageModes
		@test JEMSS.getPointsCoverageMode(demandCoverage, pcm.travelMode, pcm.coverTime) == pcm
	end
end

@testset "demand coverage calc" begin
	# compare results of fast (using sim.demandCoverage) and slow calculation of demand point coverage
	
	stationsNumAmbs = ones(Int, numStations)
	pointCoverCounts = Vector{Int}(numPoints)
	times = vcat(travel.setsStartTimes, demand.setsStartTimes) |> unique |> sort # simulating to these points in time will 
	for t in times
		# simulateToTime!(sim, t) # probably unnecessary, and sim.time may stop before t
		sim.time = t # is this safe?
		demandsPointSetsCoverCounts = JEMSS.calcPointSetsCoverCounts!(sim, sim.time, stationsNumAmbs)
		for demandPriority in priorities
			travelMode = JEMSS.getTravelMode!(travel, sim.responseTravelPriorities[demandPriority], sim.time)
			coverTime = sim.demandCoverTimes[demandPriority]
			pointCoverCounts[:] = 0
			for (i,station) in enumerate(stations)
				for (j,point) in enumerate(points)
					# calculate travel time from station to point
					time1 = offRoadTravelTime(travelMode, station.nearestNodeDist)
					pathTravelTime = shortestPathTravelTime(sim.net, travelMode.index, station.nearestNodeIndex, point.nearestNodeIndex)
					time2 = offRoadTravelTime(travelMode, point.nearestNodeDist)
					travelTime = time1 + pathTravelTime + time2
					
					# check if station covers point
					if travelTime <= coverTime
						pointCoverCounts[j] += stationsNumAmbs[i]
					end
				end
			end
			
			# compare pointSetsCoverCounts with pointCoverCounts
			pointSetsCoverCounts = demandsPointSetsCoverCounts[demandPriority] # to be checked
			pointsCoverageMode = JEMSS.getPointsCoverageMode(demandCoverage, travelMode, coverTime)
			pointSets = pointsCoverageMode.pointSets # shorthand
			@test all(i -> all(j -> pointCoverCounts[j] == pointSetsCoverCounts[i], pointSets[i]), 1:length(pointSets)) # check that calculated point cover counts are correct
			
			# also check pointSetsDemands match sum of pointDemands for points in pointSets
			demandMode = JEMSS.getDemandMode!(demand, demandPriority, sim.time)
			rasterIndex = demandMode.rasterIndex # shorthand
			pointSetsDemands = demandCoverage.pointSetsDemands[pointsCoverageMode.index, rasterIndex]
			pointDemands = demandCoverage.rastersPointDemands[rasterIndex]
			@test all(i -> isapprox(sum(j -> pointDemands[j], pointSets[i]), pointSetsDemands[i], rtol = eps(Float)*length(pointSets[i])), 1:length(pointSets)) # check that pointSetsDemands correctly stores the total demand of each pointSets
		end
	end
	
	# reset sim state, in case sim is needed again
	sim.time = sim.startTime
	travel.recentSetsStartTimesIndex = 1
	demand.recentSetsStartTimesIndex = 1
end
