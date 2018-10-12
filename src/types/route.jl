# given a route, a priority for travel, a start time, and an end location,
# (and an end node which should be the node nearest to end location),
# modify current route
function changeRoute!(sim::Simulation, route::Route, priority::Priority, startTime::Float, endLoc::Location, endFNode::Int)
	
	# shorthand:
	map = sim.map
	net = sim.net
	travel = sim.travel
	
	# get data on current route before changing
	startLoc = getRouteCurrentLocation!(net, route, startTime)
	
	travelMode = getTravelMode!(travel, priority, startTime)
	
	(startFNode, startFNodeTravelTime) = getRouteNextNode!(sim, route, travelMode.index, startTime)
	startFNodeTime = startTime + startFNodeTravelTime
	
	(pathTravelTime, rNodes) = shortestPathData(net, travelMode.index, startFNode, endFNode)
	
	# shorthand:
	fNetTravel = travelMode.fNetTravel
	fNodeFromRNodeTime = fNetTravel.fNodeFromRNodeTime
	fNodeToRNodeTime = fNetTravel.fNodeToRNodeTime
	
	## change route
	
	route.priority = priority
	route.travelModeIndex = travelMode.index
	
	# start and end fNodes and times
	route.startFNode = startFNode
	route.startFNodeTime = startFNodeTime
	route.endFNode = endFNode
	route.endFNodeTime = startFNodeTime + pathTravelTime
	
	# start and end rNodes and times
	route.startRNode = rNodes[1]
	route.endRNode = rNodes[2]
	if route.startRNode != nullIndex
		@assert(route.endRNode != nullIndex)
		route.startRNodeTime = startFNodeTime + fNodeToRNodeTime[startFNode][route.startRNode]
		route.endRNodeTime = route.endFNodeTime - fNodeFromRNodeTime[endFNode][route.endRNode]
	else
		route.startRNodeTime = nullTime
		route.endRNodeTime = nullTime
	end
	
	# start and end locations and times
	route.startLoc = startLoc
	route.startTime = startTime
	route.endLoc = endLoc
	route.endTime = route.endFNodeTime + offRoadTravelTime(travelMode, map, net.fGraph.nodes[endFNode].location, endLoc)
	
	# recent rArc, recent fNode, next fNode
	setRouteStateBeforeStartFNode!(route, startTime)
	
	# first rArc
	setRouteFirstRArc!(net, route)
end

# given a route and time, get current location
function getRouteCurrentLocation!(net::Network, route::Route, time::Float)

	updateRouteToTime!(net, route, time)
	
	fNodes = net.fGraph.nodes # shorthand
	currentLoc = nothing # init
	if time <= route.startTime
		currentLoc = route.startLoc
		
	elseif time >= route.endTime
		currentLoc = route.endLoc
		
	elseif time <= route.startFNodeTime # and time > route.startTime
		# currently between startLoc and startFNode
		currentLoc = linearInterpLocation(route.startLoc, fNodes[route.startFNode].location, route.startTime, route.startFNodeTime, time)
		
	elseif time >= route.endFNodeTime # and time < route.endTime
		# currently between endFNode and endLoc
		currentLoc = linearInterpLocation(fNodes[route.endFNode].location, route.endLoc, route.endFNodeTime, route.endTime, time)
		
	else
		# currently somewhere on network
		currentLoc = linearInterpLocation(fNodes[route.recentFNode].location, fNodes[route.nextFNode].location, route.recentFNodeTime, route.nextFNodeTime, time)
	end
	
	return currentLoc
end

# given route and current time,
# return next fNode index and time remaining to reach this node
# if already past last node in route, return that fNode and time to return to it
# testing: in certain cases, use given travel mode to determine travel time to the node
function getRouteNextNode!(sim::Simulation, route::Route, travelModeIndex::Int, time::Float)
	
	# shorthand:
	map = sim.map
	net = sim.net
	travelModes = sim.travel.modes
	
	# first need to update route
	updateRouteToTime!(net, route, time)
	
	# if route already finished, return nearest node
	# also return travel time to node based on route.endLoc
	if route.endTime <= time
		nearestFNode = route.endFNode
		travelTime = offRoadTravelTime(travelModes[travelModeIndex], map, route.endLoc, net.fGraph.nodes[nearestFNode].location)
		
		return nearestFNode, travelTime
	end
	
	nextFNode = nullIndex # init
	travelTime = nullTime # init
	if route.nextFNode != nullIndex
		# between startLoc and endFNode, go to nextFNode
		nextFNode = route.nextFNode
		travelTime = route.nextFNodeTime - time
	else
		# between endFNode and endLoc, return to recentFNode
		nextFNode = route.recentFNode
		travelTime = time - route.recentFNodeTime
	end
	
	# testing: changing travel time to match bartsim code...
	# scale travel time according to any change in travel mode
	if route.travelModeIndex != travelModeIndex
		if route.nextFNode == nullIndex
			# currently somewhere between route.endFNode and route.endLoc
			travelTime *= travelModes[route.travelModeIndex].offRoadSpeed / travelModes[travelModeIndex].offRoadSpeed
		end
	end
	
	if checkMode
		@assert(nextFNode != nullIndex)
		@assert(travelTime >= 0)
	end
	
	return nextFNode, travelTime
end

# updates route fields for given time
function updateRouteToTime!(net::Network, route::Route, time::Float)
	
	if time <= route.startFNodeTime
		@assert(route.nextFNode == route.startFNode) # see setRouteStateBeforeStartFNode!()
		return
	elseif time >= route.endFNodeTime
		setRouteStateAfterEndFNode!(route, time)
		@assert(route.recentFNode == route.endFNode) # see setRouteStateAfterEndFNode!()
		return
	else
		# currently somewhere on network
		updateRouteRecentRArc!(net, route, time)
		updateRouteRecentRArcFNode!(net, route, time)
	end
	
end

# update route.recentRArc and similar fields for given time
# assumes that route is still on network
function updateRouteRecentRArc!(net::Network, route::Route, time::Float)
	
	if route.firstRArc == nullIndex
		return # do nothing
	end
	
	# should be on network
	@assert(route.startFNodeTime <= time < route.endFNodeTime)
	
	if route.recentFNode == nullIndex
		# previously not on route, set to be at startFNode
		setRouteStateAfterStartFNode!(net, route, time)
	end
	
	if time < route.recentRArcEndTime # or: time < min(route.recentRArcEndTime, route.endFNodeTime)
		return # still on recentRArc, do nothing
	end
	
	# currently have: min(route.recentRArcEndTime, route.endFNodeTime) <= time
	# also know now that recentRArc is no longer current and so needs updating
	
	# shorthand:
	rNetTravel = net.rNetTravels[route.travelModeIndex]
	rNodeFNode = net.rNodeFNode
	
	@assert(route.startRNode != nullIndex != route.endRNode)
	
	if route.endRNodeTime <= time
		# on last rArc
		route.recentRArc = findRArcFromFNodeToFNode(net, rNodeFNode[route.endRNode], route.endFNode)
		route.recentRArcRecentFNode = 1 # remaining fNode data will be set in updateRouteRecentRArcFNode!()
		route.recentRArcStartTime = route.endRNodeTime
		route.recentRArcEndTime = route.endRNodeTime + rNetTravel.arcTimes[route.recentRArc]
		return
	end
	
	# somewhere between startRNode and endRNode
	if route.startRNode != route.endRNode
		rArcs = net.rGraph.arcs # shorthand
		rArc = rArcs[route.recentRArc]
		recentRNode = rArc.fromNodeIndex
		nextRNode = rArc.toNodeIndex
		while route.recentRArcEndTime <= time && nextRNode != route.endRNode
			recentRNode = nextRNode
			
			rArcIndex = shortestPathNextRArc(net, route.travelModeIndex, nextRNode, route.endRNode)
			nextRNode = rArcs[rArcIndex].toNodeIndex
			
			route.recentRArc = rArcIndex
			route.recentRArcStartTime = route.recentRArcEndTime
			route.recentRArcEndTime += rNetTravel.arcTimes[route.recentRArc] # do not use route.recentRArcEndTime += rNetTravel.spTimes[recentRNode, nextRNode], as this may have low precision
		end
		route.recentRArcRecentFNode = 1 # remaining fNode data will be set in updateRouteRecentRArcFNode!()
	end
	
end

# for a given route and current time,
# updates 'recentRArcRecentFNode', 'recentRArcNextFNode', and other data in route
function updateRouteRecentRArcFNode!(net::Network, route::Route, time::Float)
	# assumes that route is currently still on the recent rArc
	
	if route.firstRArc == nullIndex
		return # do nothing
	end
	
	# check that recent rArc is still current
	@assert(route.recentFNodeTime != nullTime)
	@assert(route.recentFNodeTime <= time < min(route.recentRArcEndTime, route.endFNodeTime))
	
	if time < route.nextFNodeTime
		return # do nothing
	end
	
	# shorthand:
	fNetTravel = net.fNetTravels[route.travelModeIndex]
	rArcFNodes = net.rArcFNodes[route.recentRArc]
	recentRArcFNodesTimes = fNetTravel.rArcFNodesTimes[route.recentRArc]
	
	route.recentRArcRecentFNode = findMaxIndexLeqTime(recentRArcFNodesTimes, time - route.recentRArcStartTime, route.recentRArcRecentFNode)
	@assert(route.recentRArcRecentFNode < length(rArcFNodes))
	route.recentFNode = rArcFNodes[route.recentRArcRecentFNode]
	route.recentFNodeTime = route.recentRArcStartTime + recentRArcFNodesTimes[route.recentRArcRecentFNode]
	
	route.recentRArcNextFNode = route.recentRArcRecentFNode + 1
	route.nextFNode = rArcFNodes[route.recentRArcNextFNode]
	route.nextFNodeTime = route.recentRArcStartTime + recentRArcFNodesTimes[route.recentRArcNextFNode]
end

# for a vector of strictly increasing time values,
# return the index of the largest value less than or equal (leq) to the given time
# minIndex is an optional lower bound on the index
function findMaxIndexLeqTime(times::Vector{Float}, time::Float, minIndex::Int = 1)
	
	if checkMode
		@assert(issorted(times, lt=<=)) # values should be strictly increasing
	end
	
	n = length(times)
	@assert(1 <= minIndex < n)
	if times[n] <= time
		return n
	end
	@assert(times[minIndex] <= time < times[n])
	
	# find range that contains given time
	# start at minIndex and step with increasing step sizes
	i = minIndex
	j = i + 1
	step = 2
	while times[j] <= time
		i = j
		j += step
		j >= n ? (j = n; break) : step *= 2
	end
	# should now have: times[i] <= time < times[j]
	
	# binary search between i and j
	# maintain the property: times[i] <= time < times[j]
	while i < j-1
		k = div(i+j,2)
		times[k] <= time ? i = k : j = k
	end
	
	return i
end

# given a route with the start and end fNodes and rNodes already set,
# set the firstRArc of the route
function setRouteFirstRArc!(net::Network, route::Route)
	@assert(route.startFNode != nullIndex != route.endFNode)
	
	rNodeFNode = net.rNodeFNode # shorthand
	
	if route.startFNode == route.endFNode
		# no firstRArc needed
		route.firstRArc = nullIndex
		return
		
	elseif route.startRNode == nullIndex # and route.startFNode != route.endFNode
		@assert(route.endRNode == nullIndex)
		# find single rArc which contains path from startFNode to endFNode
		route.firstRArc = findRArcFromFNodeToFNode(net, route.startFNode, route.endFNode)
		
	elseif route.startFNode != rNodeFNode[route.startRNode] # and neither nodes are null
		# find arc from startFNode to startRNode
		route.firstRArc = findRArcFromFNodeToFNode(net, route.startFNode, rNodeFNode[route.startRNode])
		
	elseif route.startRNode == route.endRNode # and route.startFNode == rNodeFNode[route.startRNode]
		@assert(rNodeFNode[route.endRNode] != route.endFNode)
		# find arc from startRNode to endFNode
		route.firstRArc = findRArcFromFNodeToFNode(net, rNodeFNode[route.startRNode], route.endFNode)
		
	else # nullIndex != route.startRNode != route.endRNode != nullIndex
		# find first arc on path from startRNode to endRNode
		route.firstRArc = shortestPathNextRArc(net, route.travelModeIndex, route.startRNode, route.endRNode)
	end
	
	@assert(route.firstRArc != nullIndex)
end

# set temporally varying fields of route to represent state before reaching startFNode
function setRouteStateBeforeStartFNode!(route::Route, time::Float)
	@assert(time <= route.startFNodeTime) # @assert(time < route.startFNodeTime) does not always work
	
	# recent rArc
	route.recentRArc = nullIndex
	route.recentRArcStartTime = nullTime
	route.recentRArcEndTime = nullTime
	
	# recent fNode
	route.recentRArcRecentFNode = nullIndex
	route.recentFNode = nullIndex
	route.recentFNodeTime = nullTime
	
	# next fNode
	route.recentRArcNextFNode = nullIndex
	route.nextFNode = route.startFNode
	route.nextFNodeTime = route.startFNodeTime
end

# set temporally varying fields of route to represent state just after leaving startFNode
function setRouteStateAfterStartFNode!(net::Network, route::Route, time::Float)
	@assert(route.startFNodeTime <= time < route.endFNodeTime)
	@assert(route.firstRArc != nullIndex)
	# or: if (route.firstRArc == nullIndex) setRouteStateBeforeStartFNode!(route, time); return; end
	
	# shorthand:
	rNetTravel = net.rNetTravels[route.travelModeIndex]
	fNetTravel = net.fNetTravels[route.travelModeIndex]
	fNodeFromRNodeTime = fNetTravel.fNodeFromRNodeTime
	fNodeToRNodeTime = fNetTravel.fNodeToRNodeTime
	rArcFNodesTimes = fNetTravel.rArcFNodesTimes
	rArcFNodes = 	net.rArcFNodes
	rArc = net.rGraph.arcs[route.firstRArc]
	firstRArcFNodesTimes = rArcFNodesTimes[route.firstRArc]
	firstRArcFNodes = rArcFNodes[route.firstRArc]
	
	# recent rArc
	route.recentRArc = route.firstRArc
	route.recentRArcStartTime = route.startFNodeTime - fNodeFromRNodeTime[route.startFNode][rArc.fromNodeIndex]
	route.recentRArcEndTime = route.recentRArcStartTime + rNetTravel.arcTimes[route.firstRArc]
	
	# recent fNode
	route.recentRArcRecentFNode = net.rArcFNodeIndex[route.firstRArc][route.startFNode]
	route.recentFNode = route.startFNode
	route.recentFNodeTime = route.startFNodeTime
	if checkMode
		@assert(firstRArcFNodes[route.recentRArcRecentFNode] == route.startFNode)
	end
	
	# next fNode
	route.recentRArcNextFNode = route.recentRArcRecentFNode + 1
	route.nextFNode = firstRArcFNodes[route.recentRArcNextFNode]
	route.nextFNodeTime = route.recentRArcStartTime + firstRArcFNodesTimes[route.recentRArcNextFNode]
end

# set temporally varying fields of route to represent state after leaving endFNode
function setRouteStateAfterEndFNode!(route::Route, time::Float)
	@assert(route.endFNodeTime <= time)
	
	# recent rArc
	route.recentRArc = nullIndex
	route.recentRArcStartTime = nullTime
	route.recentRArcEndTime = nullTime
	
	# recent fNode
	route.recentRArcRecentFNode = nullIndex
	route.recentFNode = route.endFNode
	route.recentFNodeTime = route.endFNodeTime
	
	# next fNode
	route.recentRArcNextFNode = nullIndex
	route.nextFNode = nullIndex
	route.nextFNodeTime = nullTime
end

"""
	shortestRouteTravelTime!(sim::Simulation;
		startLoc::Location, firstNode::Int, dist1::Float, time1::Float, route::Route,
		endLoc::Location, lastNode::Int, dist2::Float, time2::Float,
		travelMode::TravelMode, travelPriority::Priority, startTime::Float)
Returns the travel time of the shortest route, given information about the start and end.

# Keyword arguments
Requires one of (in order of preference):
- `route` and `startTime` (use this if considering starting on `route`)
- `firstNode` (index of first node in new route) and `time1` (travel time from `startLoc` to `firstNode`)
- `firstNode` and `dist1` (distance between `startLoc` and `firstNode`)
- `startLoc` (can additionally pass in `firstNode` or `time1`)
Also requires one of (in order of preference):
- `lastNode` (index of last node in new route) and `time2` (travel time from `lastNode` to `endLoc`)
- `lastNode` and `dist2` (distance between `lastNode` and `endLoc`)
- `endLoc` (can additionally pass in `lastNode` or `time2`)
Also requires one of (in order of preference):
- `travelMode`
- `travelPriority` and `startTime`

See also: [`shortestPathTravelTime`](@ref)
"""
# mutates: route, sim.travel
function shortestRouteTravelTime!(sim::Simulation;
	route::Union{Route, Void} = nothing, startLoc::Union{Location, Void} = nothing, firstNode::Int = nullIndex, dist1::Float = nullDist, time1::Float = nullTime,
	endLoc::Union{Location, Void} = nothing, lastNode::Int = nullIndex, dist2::Float = nullDist, time2::Float = nullTime,
	travelMode::Union{TravelMode, Void} = nothing, travelPriority::Priority = nullPriority, startTime::Float = nullTime)
	
	# determine travel mode
	if travelMode == nothing
		@assert(startTime != nullTime)
		@assert(travelPriority != nullPriority)
		travelMode = getTravelMode!(sim.travel, travelPriority, startTime)
	end
	
	# find firstNode, and time to reach it from startLoc
	fNodes = sim.net.fGraph.nodes # shorthand
	if route == nothing
		if firstNode == nullIndex
			(firstNode, dist1) = findNearestNodeInGrid(sim.map, sim.grid, fNodes, startLoc)
		end
		if time1 == nullTime
			if dist1 == nullDist
				@assert(startLoc != nothing)
				dist1 = normDist(sim.map, startLoc, fNodes[firstNode].location)
			end
			time1 = offRoadTravelTime(travelMode, dist1)
		end
	else
		@assert(startTime != nullTime)
		(firstNode, time1) = getRouteNextNode!(sim, route, travelMode.index, startTime)
	end
	
	# find lastNode, and time from it to endLoc
	if lastNode == nullIndex
		(lastNode, dist2) = findNearestNodeInGrid(sim.map, sim.grid, fNodes, endLoc)
	end
	if time2 == nullTime
		if dist2 == nullDist
			@assert(endLoc != nothing)
			dist2 = normDist(sim.map, endLoc, fNodes[lastNode].location)
		end
		time2 = offRoadTravelTime(travelMode, dist2)
	end
	
	# path (on-road) travel time
	pathTravelTime = shortestPathTravelTime(sim.net, travelMode.index, firstNode, lastNode)
	
	# total travel time
	travelTime = time1 + pathTravelTime + time2
	
	return travelTime
end
