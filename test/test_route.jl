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

runGenConfig("data/regions/small/4/gen_config.xml", overwriteOutputPath = true, doPrint = false)

@testset "route update" begin
	# Test updating of route with updateRouteToTime! function.
	# In particular, tests updating a route to a time at which it will be very close to a node.
	# Passes if assertions in routing code all pass.
	
	sim = initSim("data/regions/small/4/sim_config.xml", doPrint = false)
	
	# shorthand
	net = sim.net
	fNodes = net.fGraph.nodes
	numFNodes = length(fNodes)
	
	@assert(numFNodes > length(net.rGraph.nodes)) # check that some nodes are only in fGraph
	
	updateRoute!(route::Route, time::Float) = if (route.recentUpdateTime <= time && route.startTime <= time) JEMSS.updateRouteToTime!(net, route, time) end
	
	# test route between each pair of nodes, with and without off-road travel
	startTime = 0.0
	minArcTime = minimum(fNetTravel -> minimum(fNetTravel.arcTimes), net.fNetTravels)
	offRoadSpeed = 60*24 # 60 km/hr, though any positive value will do here
	offRoadDists = [0.0, offRoadSpeed * minArcTime]
	for i = 1:numFNodes, j = 1:numFNodes, offRoadDist in offRoadDists
		route = Route()
		offRoadTime = offRoadDist / offRoadSpeed
		JEMSS.initRoute!(route, currentLoc = fNodes[i].location, nextFNode = i, nextFNodeDist = offRoadDist)
		JEMSS.changeRoute!(sim, route, highPriority, startTime + offRoadTime, fNodes[j].location, j)
		updateRoute!(route, startTime)
		updateRoute!(route, (startTime + route.startFNodeTime) / 2)
		while route.status != routeAfterEndNode
			time = route.nextFNodeTime
			dt = (time + 1) * eps(Float)
			updateRoute!(route, time - dt) # just before next node
			updateRoute!(route, time) # just at next node
			updateRoute!(route, time + dt) # just after next node
			updateRoute!(route, time + minArcTime / 2) # at least halfway between next node future next node
		end
		updateRoute!(route, route.endTime)
	end
	@test true
end

@testset "route distance" begin
	# Test calcRouteDistance! function.
	# Compare against a calculation that relies on each arc distance being equal to the calculated distance between the arc's nodes.
	# Assumes that shortestPathDistance function is correct (this is tested in test_network.jl).
	
	sim = initSim("data/regions/small/1/sim_config.xml", doPrint = false)
	
	# shorthand
	net = sim.net
	fNodes = net.fGraph.nodes
	
	# check that arc distance is equal to distance between two nodes, measured by normDist function
	JEMSS.normDist(arc::Arc) = normDist(sim.map, fNodes[arc.fromNodeIndex].location, fNodes[arc.toNodeIndex].location)
	@assert(all(arc -> isapprox(normDist(arc), arc.distance, rtol = eps()), net.fGraph.arcs))
	
	# set absolute error tolerance between the two distance calculations
	atol = eps() * (sim.map.xRange * sim.map.xScale + sim.map.yRange * sim.map.yScale) # eps * upper bound on route distance
	
	# returns distance of route travelled up to given time,
	# assuming that arc distance is equal to calculated distance between arc's nodes
	function calcRouteDistance2(net::Network, map::Map, route::Route, time::Float)
		fNodes = net.fGraph.nodes # shorthand
		dist = 0.0
		currentLoc = JEMSS.getRouteCurrentLocation!(net, route, time)
		if route.recentFNode == nullIndex
			dist += normDist(map, route.startLoc, currentLoc) # off-road distance before startFNode
		else
			dist += normDist(map, route.startLoc, fNodes[route.startFNode].location) # off-road distance before startFNode
			dist += shortestPathDistance(net, route.travelModeIndex, route.startFNode, route.recentFNode)
			dist += normDist(map, fNodes[route.recentFNode].location, currentLoc) # off-road distance after recentFNode
		end
		return dist
	end
	
	dt = 0.0001 # time step
	endTime = min(1.0, sim.calls[end].arrivalTime)
	allPass = true
	for t = sim.startTime : dt : endTime
		simulate!(sim, time = t)
		for amb in sim.ambulances
			route = amb.route # shorthand
			if route.travelModeIndex != nullIndex
				dist1 = JEMSS.calcRouteDistance!(sim, route, t)
				dist2 = calcRouteDistance2(net, sim.map, route, t)
				allPass &= isapprox(dist1, dist2, atol = atol) # compare distances
			end
		end
	end
	@test allPass
end
