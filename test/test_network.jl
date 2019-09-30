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

using LightGraphs
using SparseArrays

function readNetworkFiles(nodesFilename::String, arcsFilename::String)
	net = Network()
	fGraph = net.fGraph # shorthand
	fGraph.nodes = readNodesFile(nodesFilename)
	(fGraph.arcs, arcTravelTimes) = readArcsFile(arcsFilename)
	
	# create map
	mp = Map()
	mp.xMin = minimum(node -> node.location.x, fGraph.nodes) - 1
	mp.xMax = maximum(node -> node.location.x, fGraph.nodes) + 1
	mp.yMin = minimum(node -> node.location.y, fGraph.nodes) - 1
	mp.yMax = maximum(node -> node.location.y, fGraph.nodes) + 1
	
	JEMSS.initGraph!(fGraph)
	JEMSS.setArcDistances!(fGraph, mp)
	JEMSS.checkGraph(fGraph, mp)
	JEMSS.initFNetTravels!(net, arcTravelTimes)
	JEMSS.createRGraphFromFGraph!(net)
	JEMSS.checkGraph(net.rGraph, mp)
	JEMSS.createRNetTravelsFromFNetTravels!(net)
	JEMSS.setCommonFNodes!(net, Int[])
	
	return net
end

function fNetTravelShortestPathData(net::Network, travelModeIndex::Int)
	# Return spTimes and spSuccs for fNetTravel:
	# spTimes[i,j] gives the travel time from fNode i to fNode j along shortest path
	# spSuccs[i,j] gives index of successor node on shortest path from fNode i to fNode j
	
	# Most of this code is copied from the function JEMSS.calcRNetTravelShortestPaths!,
	# but changed to be for fNetTravel instead of rNetTravel.
	
	# shorthand
	fGraph = net.fGraph
	numNodes = length(fGraph.nodes)
	
	# calculate all-pairs shortest-paths data for the network with Dijkstra's algorithm
	revGraph = LightGraphs.reverse(fGraph.light) # reverse graph
	revArcTimes = spzeros(FloatSpTime, numNodes, numNodes)
	for arc in fGraph.arcs
		@assert(revArcTimes[arc.toNodeIndex, arc.fromNodeIndex] == 0)
		revArcTimes[arc.toNodeIndex, arc.fromNodeIndex] = net.fNetTravels[travelModeIndex].arcTimes[arc.index]
	end
	spTimes = Array{FloatSpTime,2}(undef, numNodes, numNodes) # spTimes[i,j] is the shortest-path travel time from node i to node j, to be calculated
	spSuccs = Array{IntRNode,2}(undef, numNodes, numNodes) # spSuccs[i,j] gives the successor node of i on shortest path from node i to node j, to be calculated
	for i = 1:numNodes
		# calculate shortest paths for origin node i
		spData = LightGraphs.dijkstra_shortest_paths(revGraph, i, revArcTimes) # type: LightGraphs.DijkstraState{Float}
		spTimes[:,i] = spData.dists
		spSuccs[:,i] = spData.parents
	end
	
	# change spSuccs where needed
	# where a path i -> l has spSuccs[i,l] = k, but spSuccs[i,k] = j,
	# need to change all cases where spSuccs[i,-] = k to instead have spSuccs[i,-] = j
	# this case can happen when an arc (i,k) has the same travel time as a path (i,j,k)
	for fArc in fGraph.arcs
		i = fArc.fromNodeIndex
		j = k = fArc.toNodeIndex
		while j != spSuccs[i,j]
			j = spSuccs[i,j]
		end
		if j != k
			# for any paths from i that have successor node k, change successor to be j
			spSuccs[i, spSuccs[i,:] .== k] .= j
		end
	end
	
	return spTimes, spSuccs
end

@testset "shortest path travel time" begin
	# Check shortestPathTravelTime function.
	# Compare against a simpler (but slower) function that calculates travel times for the whole full network,
	# unlike shortestPathTravelTime which takes advantage of the reduced network.
	
	testNetworksFolder = "data/network/small"
	for testFolderName in readdir(testNetworksFolder)
		folder = joinpath(testNetworksFolder, testFolderName)
		net = readNetworkFiles(joinpath(folder, "nodes.csv"), joinpath(folder, "arcs.csv"))
		for travelModeIndex = 1:length(net.fNetTravels)
			rtol = eps(FloatSpTime)
			atol = 2 * maximum(net.rNetTravels[travelModeIndex].arcTimes) * eps(FloatSpTime) # in NetTravel type, precision of fNodeToRNodeTime and fNodeFromRNodeTime (used in sp time calculation) are affected by rArc travel time. Multiplied by 2 since this can affect the travel time at both the start and end of the path.
			
			# compare travel time from shortestPathTravelTime function with spTimes for full network
			spTimes, spSuccs = fNetTravelShortestPathData(net, travelModeIndex)
			numNodes = length(net.fGraph.nodes) # shorthand
			allPass = true
			for i = 1:numNodes, j = 1:numNodes
				t1 = shortestPathTravelTime(net, travelModeIndex, i, j)
				t2 = spTimes[i,j]
				allPass &= abs(t1-t2) <= atol + rtol*max(t1,t2) # equivalent to isapprox(t1, t2, atol = atol, rtol = rtol) for julia v0.6
			end
			@test allPass
		end
	end
end

@testset "shortest path distance" begin
	# Check shortestPathDistance function.
	# Compare against a simpler (but slower) function that calculates distances based on arc traversed in full network,
	# unlike shortestPathDistance which takes advantage of the reduced network.
	
	# Calculate distance along shortest path between two nodes, using fGraph.
	# spSuccs[i,j] gives successor of node i in shortest path from node i to node j, see fNetTravelShortestPathData.
	function shortestPathDistance2(net::Network, spSuccs::Array{IntRNode,2}, startFNode::Int, endFNode::Int)::Float
		dist = 0.0
		fNode = startFNode
		while fNode != endFNode
			nextFNode = spSuccs[fNode, endFNode]
			arcIndex = net.fGraph.nodePairArcIndex[fNode, nextFNode]
			dist += net.fGraph.arcs[arcIndex].distance
			fNode = nextFNode
		end
		return dist
	end
	
	testNetworksFolder = "data/network/small"
	for testFolderName in readdir(testNetworksFolder)
		folder = joinpath(testNetworksFolder, testFolderName)
		net = readNetworkFiles(joinpath(folder, "nodes.csv"), joinpath(folder, "arcs.csv"))
		for travelModeIndex = 1:length(net.fNetTravels)
			# compare shortestPathDistance and shortestPathDistance2 for each pair of nodes in fGraph
			spTimes, spSuccs = fNetTravelShortestPathData(net, travelModeIndex)
			numNodes = length(net.fGraph.nodes) # shorthand
			allPass = true
			for i = 1:numNodes, j = 1:numNodes
				d1 = shortestPathDistance(net, travelModeIndex, i, j) # from JEMSS
				d2 = shortestPathDistance2(net, spSuccs, i, j)
				allPass &= isapprox(d1, d2, rtol = eps())
			end
			@test allPass
		end
	end
end
