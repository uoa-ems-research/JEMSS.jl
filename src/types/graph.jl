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

# miscellaneous initialisiations for graph
function initGraph!(graph::Graph)
	# should already have nodes, arcs
	numNodes = length(graph.nodes)
	@assert(numNodes > 0 && length(graph.arcs) > 0)
	
	graph.arcDists = [arc.distance for arc in graph.arcs]
	
	graph.light = LightGraphs.DiGraph(numNodes) # create LightGraph version of graph
	for arc in graph.arcs
		LightGraphs.add_edge!(graph.light, arc.fromNodeIndex, arc.toNodeIndex) # add edges to light graph
	end
	graph.fadjList = graph.light.fadjlist
	graph.badjList = graph.light.badjlist
	
	if !graph.isReduced
		graph.nodePairArcIndex = spzeros(Int, numNodes, numNodes) # sparse, to save on memory
		for arc in graph.arcs
			@assert(graph.nodePairArcIndex[arc.fromNodeIndex, arc.toNodeIndex] == 0) # check that there is only one arc for this ordered node pair
			graph.nodePairArcIndex[arc.fromNodeIndex, arc.toNodeIndex] = arc.index
		end
	end
end

# find nearest node by evaluating distances to all nodes from given location
# should not be used except for asserting results of faster algorithms
# can set offRoadAccessRequired in order to search nodes with off-road access, or all nodes
function findNearestNode(map::Map, nodes::Vector{Node}, location::Location;
	offRoadAccessRequired::Bool = true)
	
	numNodes = length(nodes)
	
	# compare square distances
	bestSqrDist = Inf
	chosenNode = nullIndex
	sqrDist = 0.0 # init
	for i = 1:numNodes
		node = nodes[i]
		if offRoadAccessRequired && !node.offRoadAccess
			continue
		end
		sqrDist = squareDist(map, node.location, location)
		if sqrDist < bestSqrDist
			bestSqrDist = sqrDist
			chosenNode = i
		end
	end
	dist = sqrt(bestSqrDist)
	
	return chosenNode, dist
end

# Set arc distance with normDist function,
# for arcs that pass arcFilter.
function setArcDistances!(graph::Graph, map::Map;
	arcFilter::Function = (arc->isnan(arc.distance)))
	@assert(!graph.isReduced)
	for arc in filter(arcFilter, graph.arcs)
		arc.distance = normDist(map, graph.nodes[arc.fromNodeIndex].location, graph.nodes[arc.toNodeIndex].location)
	end
end

# for graph (set of nodes, arcs), check that:
# nodes are inside map borders,
# nodes are attached to at least one arc each,
# arcs are connected to two nodes,
# arcs have non-negative travel time,
# graph is strongly connected
function checkGraph(graph::Graph, map::Map)
	
	# shorthand names
	nodes = graph.nodes
	arcs = graph.arcs
	
	numNodes = length(nodes)
	numArcs = length(arcs)
	
	# check node and arc indices
	for i = 1:numNodes
		@assert(nodes[i].index == i)
	end
	for i = 1:numArcs
		@assert(arcs[i].index == i)
	end
	
	# check arc distances
	for i = 1:numArcs
		@assert(!isnan(arcs[i].distance), "arc $i has distance NaN")
		@assert(arcs[i].distance >= 0, "arc $i has negative distance")
	end
	
	# check that nodes are inside map borders
	for i = 1:numNodes
		@assert(map.xMin <= nodes[i].location.x, "node $i outside border")
		@assert(map.xMax >= nodes[i].location.x, "node $i outside border")
		@assert(map.yMin <= nodes[i].location.y, "node $i outside border")
		@assert(map.yMax >= nodes[i].location.y, "node $i outside border")
	end
	
	# check that all nodes are attached to at least one arc
	# and check that arcs are all attached to two nodes
	nodeUsed = Vector{Bool}(undef, numNodes)
	nodeUsed[1:end] .= false
	for i = 1:numArcs
		@assert(1 <= arcs[i].fromNodeIndex && arcs[i].fromNodeIndex <= numNodes, "arc $i fromNodeIndex is outside range")
		@assert(1 <= arcs[i].toNodeIndex && arcs[i].toNodeIndex <= numNodes, "arc $i toNodeIndex is outside range")
		@assert(arcs[i].fromNodeIndex != arcs[i].toNodeIndex, "arc $i fromNodeIndex and toNodeIndex are same")
		nodeUsed[arcs[i].fromNodeIndex] = true
		nodeUsed[arcs[i].toNodeIndex] = true
	end
	for i = 1:numNodes
		@assert(nodeUsed[i], "node $i not used")
	end
	
	@assert(LightGraphs.is_weakly_connected(graph.light), "graph is not weakly connected")
	@assert(LightGraphs.is_strongly_connected(graph.light), "graph is not strongly connected")
end
