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
import OpenStreetMapX
const OSM = OpenStreetMapX

include("graph_tools.jl")

# Read OSM file with road network, return nodes and arcs
# kwargs:
# - levels: highway types that will be kept, see ROAD_CLASSES in OSM for more details.
# - boundsLLA: crop map within given bounds.
# Output may need further processing, such as only keeping the largest strongly connected component, and removing duplicate arcs.
# Note that distance calculations done here are different (more accurate) than those in JEMSS.
function readOsmNetworkFile(osmFilename::String;
	levels::Union{Set{Int},Nothing} = nothing, boundsLLA::Union{OSM.Bounds{OSM.LLA},Nothing} = nothing)
	
	if levels == nothing levels = OSM.ROAD_CLASSES |> values |> Set end
	
	# read osm file
	@assert(isfile(osmFilename))
	# get osm data; this fails if the osm file does not contain bounds, e.g. <bounds minlat="1.0" minlon="2.0" maxlat="3.0" maxlon="4.0"/>
	osmData = OSM.get_map_data(osmFilename; road_levels = levels, use_cache = false) # 'use_cache = true' ignores other kwargs (e.g. road_levels)
	(nodesENU, roadways) = (osmData.nodes, osmData.roadways)
	@assert(length(nodesENU) > 0)
	@assert(length(roadways) > 0)
	@assert(typeof(first(nodesENU)[2]) == OSM.ENU)
	
	# crop map within bounds
	if boundsLLA == nothing
		boundsLLA = osmData.bounds
		@assert(boundsLLA != OSM.Bounds{OSM.LLA}(0.0, 0.0, 0.0, 0.0), "missing bounds")
	end
	nodesLLA = OSM.LLA(nodesENU, boundsLLA)
	OSM.crop!(nodesLLA, boundsLLA, roadways)
	
	# get network data from nodes and ways
	edges, class = OSM.get_edges(nodesENU, roadways) # edges = vector of tuple of vertex keys; class[i] is class (in levels) of edges[i]
	@assert(issubset(Set(class), levels))
	vertices = OSM.get_vertices(edges) # Dict{Int,Int}; vertices[k] gives index (from 1:n) of vertex key k
	weights = OSM.distance(nodesENU, edges) # weights[i] is weight of edges[i]; appears to be equivalent to Geodesy.distance, based on output
	
	# create nodes from vertices
	numNodes = length(vertices)
	nodes = Vector{Node}(undef, numNodes)
	for (key, i) in vertices
		nodes[i] = Node()
		nodes[i].index = i
		nodes[i].location.x = nodesLLA[key].lon
		nodes[i].location.y = nodesLLA[key].lat
		nodes[i].fields["osm_key"] = key
	end
	
	# create arcs from edges
	numArcs = length(edges)
	arcs = Vector{Arc}(undef, numArcs)
	for i = 1:numArcs
		v1, v2 = edges[i] # vertex keys (to, from) for edge
		arcs[i] = Arc()
		arcs[i].index = i
		arcs[i].fromNodeIndex = vertices[v1]
		arcs[i].toNodeIndex = vertices[v2]
		arcs[i].fields["osm_class"] = class[i]
		arcs[i].fields["osm_weight"] = weights[i] # edge weight = length
		arcs[i].distance = weights[i] / 1000 # convert metres to km
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
	levels::Union{Set{Int},Nothing} = nothing, classSpeeds::Vector{Dict{Int,Float}} = Dict{Int,Float}[], classOffRoadAccess::Dict{Int,Bool} = Dict{Int,Bool}(), mergeArcsFunction::Function = mergeDuplicateOsmArcs!)
	
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
	levels::Union{Set{Int},Nothing} = nothing, boundsLLA::Union{OSM.Bounds{OSM.LLA},Nothing} = nothing, classSpeeds::Vector{Dict{Int,Float}} = Dict{Int,Float}[], classOffRoadAccess::Dict{Int,Bool} = Dict{Int,Bool}(), mergeArcsFunction::Function = mergeDuplicateOsmArcs!)
	
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
