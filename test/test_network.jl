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

@testset "shortest path travel time" begin
	testNetworksFolder = "data/network/small"
	for testFolderName in readdir(testNetworksFolder)
		folder = joinpath(testNetworksFolder, testFolderName)
		net = readNetworkFiles(joinpath(folder, "nodes.csv"), joinpath(folder, "arcs.csv"))
		
		# shorthand
		fGraph = net.fGraph
		nodes = fGraph.nodes
		numNodes = length(nodes)
		fNetTravels = net.fNetTravels
		
		for travelModeIndex = 1:length(fNetTravels)
			
			# calculate all-pairs shortest-paths data for the network with Dijkstra's algorithm
			arcTimes = spzeros(FloatSpTime, numNodes, numNodes)
			for arc in fGraph.arcs
				arcTimes[arc.fromNodeIndex, arc.toNodeIndex] = fNetTravels[travelModeIndex].arcTimes[arc.index]
			end
			spTimes = Array{FloatSpTime,2}(undef, numNodes, numNodes) # spTimes[i,j] is the shortest-path travel time from node i to node j, to be calculated
			for i = 1:numNodes
				# calculate shortest paths for origin node i
				spData = LightGraphs.dijkstra_shortest_paths(fGraph.light, i, arcTimes) # type: LightGraphs.DijkstraState{Float}
				spTimes[i,:] = spData.dists
			end
			
			rtol = eps(FloatSpTime)
			atol = 2 * maximum(net.rNetTravels[travelModeIndex].arcTimes) * eps(FloatSpTime) # in NetTravel type, precision of fNodeToRNodeTime and fNodeFromRNodeTime (used in sp time calculation) are affected by rArc travel time. Multiplied by 2 since this can affect the travel time at both the start and end of the path.
			
			# compare travel time from shortestPathTravelTime function with spTimes values
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
