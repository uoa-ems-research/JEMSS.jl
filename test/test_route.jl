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
	@assert(all(arc -> isapprox(normDist(arc), arc.distance, rtol = eps(Float)), net.fGraph.arcs))
	
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
	endTime = 1.0
	allPass = true
	@assert(endTime <= sim.calls[end].arrivalTime)
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
	@test(allPass)
end
