# dmexclp - dynamic maximum expected coverage location problem

# initialise data relevant to move up
# mutates: sim.moveUpData.dmexclpData, sim.demandCoverage, sim.demandCoverTimes
function initDmexclp!(sim::Simulation;
	demandCoverTimes = Dict([p => 8/24/60 for p in instances(Priority)]), busyFraction::Float = 0.5,
	rasterCellNumPointRows::Int = 1, rasterCellNumPointCols::Int = 1)
	@assert(!sim.used)
	
	# shorthand
	numAmbs = length(sim.ambulances)
	numStations = length(sim.stations)
	
	sim.demandCoverTimes = demandCoverTimes
	dcd = sim.moveUpData.dmexclpData # shorthand
	dcd.busyFraction = busyFraction
	
	initDemandCoverage!(sim; rasterCellNumRows = rasterCellNumPointRows, rasterCellNumCols = rasterCellNumPointCols);
	
	# calculate cover benefit values, for single demand
	dcd.marginalBenefit = (busyFraction.^[0:numAmbs-1;])*(1-busyFraction)
	
	# values that will be calculated when needed
	dcd.stationNumIdleAmbs = zeros(Int, numStations)
	dcd.stationMarginalCoverages = zeros(Float, numStations) # stationMarginalCoverages[i] gives extra coverage provided from placing newly idle ambulance at station i
	# dcd.pointSetsCoverCounts = [zeros(Int, length(pointsCoverageMode.pointSets)) for pointsCoverageMode in sim.demandCoverage.pointsCoverageModes] # pointSetsCoverCounts[i][j] = number of idle ambulances covering node set j, for demand.pointsCoverageModes i
end

function dmexclpMoveUp(sim::Simulation, newlyIdleAmb::Ambulance)
	@assert(sim.moveUpData.useMoveUp)
	@assert(newlyIdleAmb.status != ambGoingToCall)
	
	# shorthand names:
	dcd = sim.moveUpData.dmexclpData
	ambulances = sim.ambulances
	stations = sim.stations
	numStations = length(stations)
	
	# calculate the number of idle ambulances at (or travelling to) each station
	dcd.stationNumIdleAmbs[:] = 0
	for (i,ambulance) in enumerate(ambulances)
		# do not count newly idle ambulance, it has not been assigned a station
		if isAmbAvailableForMoveUp(ambulance) && i != newlyIdleAmb.index
			dcd.stationNumIdleAmbs[ambulance.stationIndex] += 1
		end
	end
	
	# ignoring newly idle amb, count number of ambulances covering each point set
	demandsPointSetsCoverCounts = calcPointSetsCoverCounts!(sim, sim.time, dcd.stationNumIdleAmbs)
	
	# find station allocation for newly idle ambulance that gives greatest
	# increase in expected demand coverage
	dcd.stationMarginalCoverages[:] = 0.0
	for demandPriority in setdiff([instances(Priority)...], [nullPriority])
		pointSetsCoverCounts = demandsPointSetsCoverCounts[demandPriority]
		pointsCoverageMode = getPointsCoverageMode!(sim, demandPriority, sim.time)
		demandMode = getDemandMode!(sim.demand, demandPriority, sim.time)
		pointSetsDemands = sim.demandCoverage.pointSetsDemands[pointsCoverageMode.index, demandMode.rasterIndex]
		# to do: if marginal coverage has already been calculated for the combination of pointsCoverageMode and rasterIndex (but for a different demand priority), the marginal coverage value should be reused (accounting for difference in demandMode.rasterMultiplier values) to save on computation.
		# to do: allow for each demand priority to have a different weight in the coverage calculation
		for i = 1:length(pointsCoverageMode.pointSets)
			pointSetMarginalCoverage = pointSetsDemands[i] * dcd.marginalBenefit[pointSetsCoverCounts[i]+1] * demandMode.rasterMultiplier # this puts equal weight on each demand priority
			for j in pointsCoverageMode.stationSets[i]
				dcd.stationMarginalCoverages[j] += pointSetMarginalCoverage
			end
		end
	end
	(bestMarginalCoverage, bestStationIndex) = findmax(dcd.stationMarginalCoverages)
	
	return [newlyIdleAmb], [stations[bestStationIndex]]
end
