# dmexclp - dynamic maximum expected coverage location problem

# initialise data relevant to move up
function initDmexclp!(sim::Simulation;
	coverTime::Float = 8/(24*60), busyFraction::Float = 0.5, demandRasterFilename::String = "")
	
	assert(isfile(demandRasterFilename))
	
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
	
	warn("Have not yet accounted for travel time from each demand location to the nearest node.")
	
	# read demand raster, set node demands
	raster = dcd.demandRaster = readRasterFile(demandRasterFilename) # shorthand
	dcd.nodeDemand = zeros(Float, numNodes)
	for i = 1:raster.nx, j = 1:raster.ny
		if raster.z[i,j] != 0
			location = Location(raster.x[i], raster.y[j])
			(nearestNodeIndex, dist) = findNearestNodeInGrid(sim.map, sim.grid, fGraph.nodes, location)
			dcd.nodeDemand[nearestNodeIndex] += raster.z[i,j]
		end
	end
	
	# determine node coverage
	dcd.stationCoverNode = Array{Bool,2}(numStations, numNodes)
	dcd.stationCoverNode[:] = false
	assert(travel.numSets == 1) # otherwise, would need coverage and busy fraction to change with time
	travelMode = getTravelMode!(travel, highPriority, sim.time) # or lowPriority?
	for i = 1:numStations
		station = sim.stations[i] # shorthand
		
		# get travel time from station to nearest node
		(node1, dist1) = (station.nearestNodeIndex, station.nearestNodeDist)
		time1 = offRoadTravelTime(travelMode, dist1) # time to reach nearest node
		
		for j = 1:numNodes
			# get shortest path travel time
			(travelTime, rNodes) = shortestPathTravelTime(net, travelMode.index, node1, j)
			if time1 + travelTime <= dcd.coverTime
				dcd.stationCoverNode[i,j] = true
			end
		end
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
