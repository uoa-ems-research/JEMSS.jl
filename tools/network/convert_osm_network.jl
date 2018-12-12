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

# for converting OpenStreetMap files with road networks to format usable for JEMSS

using JEMSS
using Geodesy # need v0.0.1 for type Bounds
import OpenStreetMap
const OSM = OpenStreetMap

include("graph_tools.jl")

# Read OSM file with road network, return nodes and arcs
# kwargs:
# - levels: highway types that will be kept, see ROAD_CLASSES in OSM for more details.
# - boundsLLA: crop map within given bounds.
# Output may need further processing, such as only keeping the largest strongly connected component, and removing duplicate arcs.
# Note that distance calculations done here are different (more accurate) than those in JEMSS.
function readOsmNetworkFile(osmFilename::String;
	levels::Union{Set{Int},Void} = nothing, boundsLLA::Union{Geodesy.Bounds{LLA},Void} = nothing)
	
	# read osm file
	@assert(isfile(osmFilename))
	(nodesLLA, highways, buildings, features) = OSM.getOSMData(osmFilename)
	
	if boundsLLA == nothing
		try
			# get map bounds from file (may be missing bounds)
			xdoc = OSM.parseMapXML(osmFilename)
			boundsLLA = OSM.getBounds(xdoc)
		catch
			# set map bounds to fit nodesLLA
			nodeCoords = collect(values(nodesLLA))
			minLat = minimum(c->c.lat, nodeCoords)
			maxLat = maximum(c->c.lat, nodeCoords)
			minLon = minimum(c->c.lon, nodeCoords)
			maxLon = maximum(c->c.lon, nodeCoords)
			boundsLLA = Geodesy.Bounds{LLA}(minLat, maxLat, minLon, maxLon)
		end
	end
	# crop map within bounds
	# even if cropping is not required, it is useful to remove highway references to nodes that may not exist
	OSM.cropMap!(nodesLLA, boundsLLA; highways=highways, buildings=buildings, features=features)
	
	# convert from LLA to ENU coordinates
	lla_reference = Geodesy.center(boundsLLA)
	nodesENU = ENU(nodesLLA, lla_reference)
	
	# create network
	if levels == nothing levels = OSM.ROAD_CLASSES |> values |> Set end
	classes = OSM.roadways(highways)
	network = OSM.createGraph(nodesENU, highways, classes, levels)
	
	graph = network.g # graph type: Graphs.GenericIncidenceList
	@assert(graph.is_directed)
	edges = OSM.getEdges(network)
	
	# create nodes from graph.vertices
	numNodes = length(graph.vertices)
	nodes = Vector{Node}(numNodes)
	for i = 1:numNodes
		vertex = graph.vertices[i]
		nodes[i] = Node()
		nodes[i].index = i
		nodes[i].location.x = nodesLLA[vertex.key].lon
		nodes[i].location.y = nodesLLA[vertex.key].lat
		nodes[i].fields["osm_key"] = vertex.key
	end
	
	# create arcs from edges
	numArcs = graph.nedges
	arcs = Vector{Arc}(numArcs)
	for i = 1:numArcs
		edge = edges[i]
		arcs[i] = Arc()
		arcs[i].index = i
		arcs[i].fromNodeIndex = edge.source.index
		arcs[i].toNodeIndex = edge.target.index
		arcs[i].fields["osm_class"] = network.class[edge.index]
		arcs[i].fields["osm_weight"] = network.w[edge.index] # edge weight = length
		@assert(network.w[edge.index] == distance(nodesENU, edge.source.key, edge.target.key))
	end
	
	return nodes, arcs
end

# merge arc2 into arc1, keeping minimum class (should have higher speed) in arc1
# arc2 will need to be removed in separate step
function mergeDuplicateOsmArcs!(arc1::Arc, arc2::Arc; classFieldName::String = "osm_class")
	@assert(arc1.fromNodeIndex == arc2.fromNodeIndex)
	@assert(arc1.toNodeIndex == arc2.toNodeIndex)
	arc1.fields[classFieldName] = min(arc1.fields[classFieldName], arc2.fields[classFieldName])
	# arcs should already have same weight (length)
end
mergeDuplicateOsmArcs!(arcs::Vector{Arc}, i::Int, j::Int) = mergeDuplicateOsmArcs!(arcs[i], arcs[j])

# Process nodes and arcs from OSM, processing performed:
# - only keep specified highway types (levels)
# - keep only the largest strongly connected component
# - merge duplicate arcs (same start and end node) according to kwarg mergeArcsFunction
# - set offRoadAccess field of nodes, according to kwarg classOffRoadAccess
# - calculate travel time for each arc
# kwargs:
# - levels: highway types that will be kept, see ROAD_CLASSES in OSM for more details.
# - classSpeeds: classSpeeds[i] has dict of the speed (km/hr) of each class of arc, for travel mode i
# - classOffRoadAccess: dict indicating which nodes of different arc classes can be used to get on and off road; node will not have off-road access if all connecting arcs have class returning false in dict
# - mergeArcsFunction: function to merge any duplicate arcs found, given the arcs and two arc indices
function convertOsmNetwork!(nodes::Vector{Node}, arcs::Vector{Arc};
	levels::Union{Set{Int},Void} = nothing, classSpeeds::Vector{Dict{Int,Float}} = [], classOffRoadAccess::Dict{Int,Bool} = Dict{Int,Bool}(), mergeArcsFunction::Function = mergeDuplicateOsmArcs!)
	
	numTravelModes = length(classSpeeds)
	
	if levels == nothing levels = OSM.ROAD_CLASSES |> values |> Set end
	
	# if classOffRoadAccess not given, assume it is true for all levels
	if classOffRoadAccess == Dict()
		classOffRoadAccess = Dict([level => true for level in levels])
	end
	
	# check that input dicts contain all levels
	for level in levels
		@assert(all(c -> haskey(c, level), classSpeeds))
		@assert(haskey(classOffRoadAccess, level))
	end
	
	# convert classSpeeds (km/hr) to classInvSpeeds (days/metre)
	# this is needed as the arc lengths (arc field "osm_weight") are in metres
	classInvSpeeds = deepcopy(classSpeeds) # will be days / metre
	for i = 1:numTravelModes
		for (class, speed) in classSpeeds[i]
			classInvSpeeds[i][class] = 1 / (speed * 1000 * 24) # days / metre
		end
	end
	
	# edit nodes and arcs as needed
	graphRemoveElts!(nodes, arcs, arcFilter = arc -> in(arc.fields["osm_class"], levels)) # keep only specified levels
	graphRemoveDisconnectedArcs!(nodes, arcs) # probably not needed, but use just in case
	graphKeepLargestComponent!(nodes, arcs) # keep only the largest strongly connected component
	graphMergeDuplicateArcs!(nodes, arcs; mergeArcsFunction = mergeArcsFunction)
	
	# set node offRoadAccess values
	for node in nodes
		node.offRoadAccess = false
	end
	for arc in arcs
		if classOffRoadAccess[arc.fields["osm_class"]]
			nodes[arc.fromNodeIndex].offRoadAccess = true
			nodes[arc.toNodeIndex].offRoadAccess = true
		end
	end
	
	# set arc travel times
	for i = 1:numTravelModes
		modeFieldName = travelModeString(i)
		for arc in arcs
			class = arc.fields["osm_class"]
			arc.fields[modeFieldName] = arc.fields["osm_weight"] * classInvSpeeds[i][class]
		end
	end
end

function convertOsmNetworkExample!(nodes::Vector{Node}, arcs::Vector{Arc})
	# the 'levels' of highways that will be kept, see ROAD_CLASSES in OSM for more details
	# e.g., levels = Set(1:5) are motorway to tertiary, including links
	# levels = Set(1:6) # roads 1:6 are motorway to residential and unclassified, including links
	levels = Set(1:6)
	
	# arc speed (km/hr) for each travel mode and each arc class
	v = [
		[110, 110, 60, 60, 60, 60, 60, 60], # travel mode 1, classes 1:8
		[ 90, 90, 45, 45, 45, 45, 45, 45] # travel mode 2, classes 1:8
	]
	classSpeeds = [Dict([j => Float(v[i][j]) for j = 1:length(v[i])]) for i = 1:length(v)]
	
	# whether nodes along road can be used to get on-road and off-road
	# where set to true, will set both ends of arc to be accessible
	classOffRoadAccess = Dict(
		1 => false, # motorway + link
		2 => false, # trunk + link
		3 => true, # primary + link
		4 => true, # secondary + link
		5 => true, # tertiary + link
		6 => true, # residential + link
	)
	
	convertOsmNetwork!(nodes, arcs; levels=levels, classSpeeds=classSpeeds, classOffRoadAccess=classOffRoadAccess)
end

function convertOsmNetworkFile(osmFilename::String;
	levels::Union{Set{Int},Void} = nothing, boundsLLA::Union{Geodesy.Bounds{LLA},Void} = nothing, classSpeeds::Vector{Dict{Int,Float}} = [], classOffRoadAccess::Dict{Int,Bool} = Dict{Int,Bool}(), mergeArcsFunction::Function = mergeDuplicateOsmArcs!)
	
	(nodes, arcs) = readOsmNetworkFile(osmFilename; levels = levels, boundsLLA = boundsLLA)
	convertOsmNetwork!(nodes, arcs; levels = levels, classSpeeds = classSpeeds, classOffRoadAccess = classOffRoadAccess, mergeArcsFunction = mergeArcsFunction)
	return nodes, arcs
end

# example use of convertOsmNetwork function
function convertOsmNetworkFileExample(osmFilename::String)
	levels = Set(1:6) # roads 1:6 are motorway to residential and unclassified, including links
	(nodes, arcs) = readOsmNetworkFile(osmFilename; levels = levels)
	convertOsmNetworkExample!(nodes, arcs)
	return nodes, arcs
end
