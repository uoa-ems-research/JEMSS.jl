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

@info("Todo: use distance calculated in changeRoute! function.")
@info("Todo: account for distance of last route at end of sim (currently not calculated).")
@info("Todo: change calcRouteDistance! function to work before ambulance route is first changed with changeRoute! function.")

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
	startFNodeDist = route.nextFNodeDist = getRouteNextNodeDist!(sim, route, startTime)
	
	(pathTravelTime, rNodes) = shortestPathData(net, travelMode.index, startFNode, endFNode)
	
	dist = route.travelModeIndex == nullIndex ? 0.0 : calcRouteDistance!(sim, route, startTime) # calcRouteDistance! should be called before resetting route.nextFNodeDist to nullDist and before setting route.startFNodeDist for new route
	
	# shorthand:
	fNetTravel = travelMode.fNetTravel
	fNodeFromRNodeTime = fNetTravel.fNodeFromRNodeTime
	fNodeToRNodeTime = fNetTravel.fNodeToRNodeTime
	
	## change route
	
	route.priority = priority
	route.travelModeIndex = travelMode.index
	
	# start and end fNodes, times, and distances
	route.startFNode = startFNode
	route.startFNodeTime = startFNodeTime
	route.startFNodeDist = startFNodeDist
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
	
	# recent rArc, recent fNode, next fNode, status
	setRouteStateBeforeStartFNode!(route, startTime)
	route.nextFNodeDist = nullDist
	
	# first rArc
	setRouteFirstRArc!(net, route)
end

# Initialise an empty route to be at a given location,
# along with the start node to travel to and distance to that node.
function initRoute!(sim::Simulation, route::Route;
	startLoc::Location = Location(), startFNode::Int = nullIndex, startFNodeDist::Float = nullDist)
	
	@assert(route.status == routeNullStatus)
	@assert(!isSameLocation(startLoc, Location()))
	@assert(startFNode != nullIndex)
	@assert(startFNodeDist >= 0.0)
	
	# make route that starts at time = Inf
	route.startTime = Inf
	route.startLoc = startLoc
	route.startFNode = startFNode
	route.startFNodeDist = startFNodeDist
	route.startFNodeTime = Inf
	# route.endTime = Inf # leave as nullTime, for getRouteNextNode!
	route.endLoc = startLoc # needed for animation
	route.endFNode = startFNode
	route.endFNodeTime = Inf
	route.nextFNode = startFNode
	setRouteStateBeforeStartFNode!(route, Inf)
	@assert(route.status == routeBeforeStartNode)
	
	# check that route functions return expected values
	if checkMode
		t = route.recentUpdateTime = 0.0
		travelModeIndex = 1
		@assert(isRouteUpToDate(route, t))
		@assert(isSameLocation(getRouteCurrentLocation!(sim.net, route, t), route.startLoc))
		@assert(getRouteNextNode!(sim, route, travelModeIndex, t)[1] == startFNode)
		@assert(getRouteNextNodeDist!(sim, route, t) == startFNodeDist)
		@assert(calcRouteDistance!(sim, route, t) == 0)
		route.recentUpdateTime = nullTime # reset
	end
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

# Get the distance of the route to the next node at the given time.
# If already past last node in route, return the distance to return to it.
# See also: getRouteNextNode!
function getRouteNextNodeDist!(sim::Simulation, route::Route, time::Float)
	# shorthand:
	map = sim.map
	net = sim.net
	fGraph = net.fGraph
	
	# first need to update route
	updateRouteToTime!(net, route, time)
	
	dist = nullDist # init
	if route.status == routeAfterEndNode
		# off-road, after last node
		dist = normDist(map, route.endLoc, fGraph.nodes[route.endFNode].location)
		if time < route.endTime
			dist *= (time - route.endFNodeTime) / (route.endTime - route.endFNodeTime)
		end
	elseif route.status == routeBeforeStartNode
		# before first node, but may be off-road or partway along first arc
		dist = route.startFNodeDist # can be different from normDist between startLoc and startFNode if starting along arc
		if time > route.startTime
			dist *= (route.startFNodeTime - time) / (route.startFNodeTime - route.startTime)
		end
	else
		@assert(route.status == routeOnPath)
		arcIndex = fGraph.nodePairArcIndex[route.recentFNode, route.nextFNode]
		@assert(arcIndex != 0)
		arc = fGraph.arcs[arcIndex]
		dist = arc.distance * (route.nextFNodeTime - time) / (route.nextFNodeTime - route.recentFNodeTime)
	end
	@assert(dist >= 0)
	
	return dist
end

function isRouteUpToDate(route::Route, time::Float)
	@assert(time != nullTime)
	@assert(route.recentUpdateTime <= time) # route should not be ahead of time
	
	# return time == route.recentUpdateTime # basic check, assuming that other fields of route are correct
	if time != route.recentUpdateTime return false end
	
	# further checks
	result = true
	if time >= route.endFNodeTime
		result &= (route.status == routeAfterEndNode)
		result &= (route.recentFNode == route.endFNode && route.nextFNode == nullIndex)
	elseif time <= route.startFNodeTime
		result &= (route.status == routeBeforeStartNode)
		result &= (route.recentFNode == nullIndex && route.nextFNode == route.startFNode)
	else # route.startFNodeTime <= time < route.endFNodeTime
		result &= (route.status == routeOnPath)
		result &= (route.recentFNode != nullIndex && route.nextFNode != nullIndex)
	end
	return result
end

# updates route fields for given time
function updateRouteToTime!(net::Network, route::Route, time::Float)
	@assert(time != nullTime)
	@assert(route.recentUpdateTime <= time)
	
	if isRouteUpToDate(route, time) return end # already up to date
	route.recentUpdateTime = time
	
	if time <= route.startFNodeTime # time < route.startFNodeTime does not always work; see setRouteStateBeforeStartFNode!()
		@assert(route.status == routeBeforeStartNode && route.nextFNode == route.startFNode) # see setRouteStateBeforeStartFNode!()
	elseif time >= route.endFNodeTime
		setRouteStateAfterEndFNode!(route, time)
		# @assert(route.status == routeAfterEndNode && route.recentFNode == route.endFNode) # see setRouteStateAfterEndFNode!()
	else
		# currently somewhere on network
		updateRouteRecentRArc!(net, route, time)
		updateRouteRecentRArcFNode!(net, route, time)
		@assert(route.status == routeOnPath)
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
	
	if route.status != routeOnPath
		# previously not on path, set to be at startFNode
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
	
	if route.endRNodeTime <= time || isapprox(route.endRNodeTime, time) # allow for imprecision in route.endRNodeTime
		# on last rArc
		route.recentRArc = findRArcFromFNodeToFNode(net, rNodeFNode[route.endRNode], route.endFNode)
		route.recentRArcRecentFNode = 1 # remaining fNode data will be set in updateRouteRecentRArcFNode!()
		route.endRNodeTime = min(route.endRNodeTime, time) # ensure that route.endRNodeTime <= time
		route.recentRArcStartTime = route.endRNodeTime
		route.recentRArcEndTime = route.endRNodeTime + rNetTravel.arcTimes[route.recentRArc]
		@assert(route.endFNodeTime < route.recentRArcEndTime) # should leave network before rArc ends
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
		if rNodeFNode[nextRNode] == route.endFNode
			@assert(isapprox(route.recentRArcEndTime, route.endFNodeTime))
			route.recentRArcEndTime = route.endFNodeTime # adjust value, to compensate for problems with numerical precision
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
	if route.nextFNode == route.endFNode
		@assert(isapprox(route.nextFNodeTime, route.endFNodeTime))
		route.nextFNodeTime = route.endFNodeTime # adjust value, to compensate for problems with numerical precision
	end
end

# for a vector of strictly increasing time values,
# return the index of the largest value less than or equal (leq) to the given time
# minIndex is an optional lower bound on the index
function findMaxIndexLeqTime(times::Vector{Float}, time::Float, minIndex::Int = 1)
	
	if checkMode && false # skip this check, it is slow
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
	
	# # or, instead of finding an upper bound for j, just set j to n
	# i = minIndex
	# j = n
	
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
	@assert(time <= route.startFNodeTime) # @assert(time < route.startFNodeTime) does not work if starting directly at route.startFNode
	
	route.status = routeBeforeStartNode
	
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
	rArc = net.rGraph.arcs[route.firstRArc]
	firstRArcFNodesTimes = fNetTravel.rArcFNodesTimes[route.firstRArc]
	firstRArcFNodes = net.rArcFNodes[route.firstRArc]
	
	route.status = routeOnPath
	
	# recent rArc
	route.recentRArc = route.firstRArc
	route.recentRArcStartTime = route.startFNodeTime - fNetTravel.fNodeFromRNodeTime[route.startFNode][rArc.fromNodeIndex]
	route.recentRArcEndTime = route.recentRArcStartTime + rNetTravel.arcTimes[route.firstRArc]
	if firstRArcFNodes[end] == route.endFNode
		@assert(isapprox(route.recentRArcEndTime, route.endFNodeTime))
		route.recentRArcEndTime = route.endFNodeTime # adjust value, to compensate for problems with numerical precision
	end
	
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
	if route.nextFNode == route.endFNode
		@assert(isapprox(route.nextFNodeTime, route.endFNodeTime))
		route.nextFNodeTime = route.endFNodeTime # adjust value, to compensate for problems with numerical precision
	end
end

# set temporally varying fields of route to represent state after leaving endFNode
function setRouteStateAfterEndFNode!(route::Route, time::Float)
	@assert(route.endFNodeTime <= time)
	
	route.status = routeAfterEndNode
	
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

Mutates: `sim.travel`, `route`

See also: [`shortestPathTravelTime`](@ref)
"""
function shortestRouteTravelTime!(sim::Simulation;
	route::Union{Route, Nothing} = nothing, startLoc::Union{Location, Nothing} = nothing, firstNode::Int = nullIndex, dist1::Float = nullDist, time1::Float = nullTime,
	endLoc::Union{Location, Nothing} = nothing, lastNode::Int = nullIndex, dist2::Float = nullDist, time2::Float = nullTime,
	travelMode::Union{TravelMode, Nothing} = nothing, travelPriority::Priority = nullPriority, startTime::Float = nullTime)
	
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

# Return the distance travelled along a route,
# from the route start time to the given time.
function calcRouteDistance!(sim::Simulation, route::Route, time::Float)::Float
	@assert(route.travelModeIndex != nullIndex)
	@assert(route.status != routeNullStatus)
	
	# shorthand
	net = sim.net
	fGraph = net.fGraph
	
	# first need to update route
	updateRouteToTime!(net, route, time)
	
	dist = 0.0
	nextFNodeDist = route.nextFNodeDist != nullDist ? route.nextFNodeDist : getRouteNextNodeDist!(sim, route, time)
	if route.status == routeBeforeStartNode
		dist = route.startFNodeDist - nextFNodeDist
	elseif route.status == routeAfterEndNode
		dist += route.startFNodeDist
		dist += shortestPathDistance(net, route.travelModeIndex, route.startFNode, route.endFNode, startRNode = route.startRNode, endRNode = route.endRNode)
		dist += nextFNodeDist # distance to return back to route.endFNode
	else
		@assert(route.status == routeOnPath)
		dist += route.startFNodeDist
		dist += shortestPathDistance(net, route.travelModeIndex, route.startFNode, route.recentFNode)
		fArcIndex = fGraph.nodePairArcIndex[route.recentFNode, route.nextFNode]
		fArc = fGraph.arcs[fArcIndex]
		dist += fArc.distance - nextFNodeDist
	end
	@assert(dist >= 0)
	
	return dist
end
