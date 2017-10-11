# dmexclp - dynamic maximum expected coverage location problem

# initialise data relevant to move up
function initDmexclp!(sim::Simulation;
	coverTime::Float = 8/(24*60), busyFraction::Float = 0.5)
	
	# shorthand names:
	dcd = sim.moveUpData.dmexclpData
	net = sim.net
	travel = sim.travel
	fGraph = net.fGraph
	
	dcd.coverTime = coverTime # (days)
	dcd.busyFraction = busyFraction
	
	numAmbs = length(sim.ambulances)
	numStations = length(sim.stations)
	numNodes = length(fGraph.nodes)
	numCalls = length(sim.calls)
	
	# calculate cover benefit values, for single demand
	q = dcd.busyFraction # shorthand
	dcd.marginalBenefit = (q.^[0:numAmbs-1;])*(1-q)
	
	# function requires some changes:
	# - should require input of demand locations and demands, rather than using actual future demand
	# - need to use network to calculate travel times to demand locations
	error("dmexclp not complete")
	
	# testing: get demand at each node by adding number of future calls
	dcd.nodeDemand = Vector{Int}(numNodes)
	dcd.nodeDemand[:] = 0
	for i = 1:numCalls
		dcd.nodeDemand[sim.calls[i].nearestNodeIndex] += 1
	end
	
	# testing: determine node coverage
	dcd.stationCoverNode = Array{Bool,2}(numStations, numNodes)
	dcd.stationCoverNode[:] = false
	assert(travel.numSets == 1) # otherwise, would need coverage and busy fraction to change with time
	travelMode = getTravelMode!(travel, lowPriority, sim.time)
	# spTimes = travelMode.rNetTravel.spTimes
	for i = 1:numStations
		startNode = sim.stations[i].nearestNodeIndex
		startTime = offRoadTravelTime(travelMode, sim.stations[i].nearestNodeDist) # time to reach node nearest station
		
		dcd.stationCoverNode[i,:] = ((spTimes[startNode,:] + startTime) .<= dcd.coverTime)
	end
	dcd.stationCoverNodes = Vector{Vector{Int}}(numStations)
	for i = 1:numStations
		dcd.stationCoverNodes[i] = find(dcd.stationCoverNode[i,:])
	end
	
	# testing: will set this to 0 for now, and calculate it each time it is needed
	dcd.stationNumIdleAmbs = Vector{Int}(numStations)
	dcd.stationNumIdleAmbs[:] = 0
	
	dcd.nodeCoverCount = Vector{Int}(numNodes)
	dcd.nodeCoverCount[:] = 0
end

function dmexclpMoveUp(sim::Simulation, newlyIdleAmb::Ambulance)
	assert(sim.moveUpData.useMoveUp)
	assert(newlyIdleAmb.status != ambGoingToCall)
	
	# shorthand names:
	ambulances = sim.ambulances
	stations = sim.stations
	dcd = sim.moveUpData.dmexclpData
	
	numAmbs = length(ambulances)
	numStations = length(stations)
	
	# calculate the number of idle ambulances at (or travelling to) each station
	dcd.stationNumIdleAmbs[:] = 0
	for i = 1:numAmbs
		# do not count newly idle ambulance, it has not been assigned a station
		if isAmbAvailableForMoveUp(ambulances[i]) && i != newlyIdleAmb.index
			dcd.stationNumIdleAmbs[ambulances[i].stationIndex] += 1
		end
	end
	
	# ignoring new idle amb, count number of ambulances covering each node
	dcd.nodeCoverCount[:] = 0
	for i = 1:numStations
		for j = dcd.stationCoverNodes[i]
			dcd.nodeCoverCount[j] += dcd.stationNumIdleAmbs[i]
		end
	end
	
	# find station allocation for new idle ambulance that gives greatest
	# increase in expected demand coverage
	extraCover = 0 # init
	bestExtraCover = -1
	bestStationIndex = nullIndex # init
	for i = 1:numStations
		# calculate additional expected coverage if new idle amb placed at station i
		extraCover = 0
		for j = dcd.stationCoverNodes[i]
			if dcd.nodeDemand[j] > 0
				extraCover += dcd.nodeDemand[j] * dcd.marginalBenefit[dcd.nodeCoverCount[j]+1]
			end
		end
		# slow:
		# for j = dcd.stationCoverNodes[i]
			# extraCover += dcd.nodeDemand[j] * dcd.marginalBenefit[dcd.nodeCoverCount[j]+1]
		# end
		
		if extraCover > bestExtraCover
			bestExtraCover = extraCover
			bestStationIndex = i
		end
	end
	
	return [newlyIdleAmb], [stations[bestStationIndex]]
end
