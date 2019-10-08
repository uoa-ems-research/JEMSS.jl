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

using Parameters

# initialise sim
testRegionDataFolder = "data/regions/small/1"
global simDemandCoverage = initSim(joinpath(testRegionDataFolder, "sim_config.xml"));

# initialise demand and demand coverage
initDemand!(simDemandCoverage, demandFilename = joinpath(testRegionDataFolder, "demand", "demand.csv"))
coverTimes = Dict([p => 20/60/24 for p in priorities])
coverTimes[highPriority] = 12/60/24
initDemandCoverage!(simDemandCoverage, coverTimes = coverTimes, rasterCellNumRows = 2, rasterCellNumCols = 2)

@testset "demand coverage init" begin
	# check sim.demandCoverage after being initialised with function initDemandCoverage!
	
	# shorthand
	global simDemandCoverage
	sim = simDemandCoverage
	@unpack demand, demandCoverage, stations, numStations, travel = sim;
	@unpack nodesPoints, points, pointsCoverageModes, pointsCoverageModeLookup = demandCoverage;
	numPoints = length(points)
	numNodes = length(sim.net.fGraph.nodes)
	
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
		@unpack pointSets, stationSets, stationsCoverPointSets, travelMode = pointsCoverageMode
		
		@assert(length(pointSets) == length(stationSets))
		@test vcat(pointSets...) |> sort == [1:numPoints;] # each point should be included in exactly one point set
		@test allunique(stationSets) # each set of stations should be unique, otherwise any duplicates should be combined
		
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
		
		@test all(i -> all(j -> findall(stationsCoverPoints[:,j]) == stationSets[i], pointSets[i]), 1:length(pointSets)) # assume that both vectors (stationsCoverPoints[:,j] and stationSets[i]) are sorted the same
		@test all(i -> Set(vcat(pointSets[stationsCoverPointSets[i]]...)) == Set(findall(stationsCoverPoints[i,:])), 1:numStations) # check stationsCoverPointSets
	end
	
	# check pointsCoverageModeLookup
	for pcm in pointsCoverageModes
		@test JEMSS.getPointsCoverageMode(demandCoverage, pcm.travelMode, pcm.coverTime) == pcm
	end
end

@testset "demand coverage calc" begin
	# compare results of fast (using sim.demandCoverage) and slow calculation of demand point coverage
	
	# shorthand
	global simDemandCoverage
	sim = simDemandCoverage
	@unpack demand, demandCoverage, stations, numStations, travel = sim;
	@unpack points = demandCoverage;
	numPoints = length(points)
	
	stationsNumAmbs = ones(Int, numStations)
	pointCoverCounts = Vector{Int}(undef, numPoints)
	times = vcat(travel.setsStartTimes, demand.setsStartTimes) |> unique |> sort # simulating to these points in time will cause each combination of travel and demand states to be checked
	for t in times
		# simulateToTime!(sim, t) # probably unnecessary, and sim.time may stop before t
		sim.time = t # is this safe?
		demandsPointSetsCoverCounts = JEMSS.calcPointSetsCoverCounts!(sim, sim.time, stationsNumAmbs)
		for demandPriority in priorities
			travelMode = JEMSS.getTravelMode!(travel, sim.responseTravelPriorities[demandPriority], sim.time)
			coverTime = demandCoverage.coverTimes[demandPriority]
			pointCoverCounts[:] .= 0
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
