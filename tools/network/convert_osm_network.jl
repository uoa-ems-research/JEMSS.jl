# for converting OpenStreetMap files with road networks to format usable for JEMSS

using JEMSS
import OpenStreetMap
using Geodesy

const OSM = OpenStreetMap

# read osm file with road network, convert to format for use by JEMSS
# processing performed:
# - only keep specified highway types (levels)
# - keep only the largest strongly connected component
# - remove duplicate arcs (same start and end node), keep arcs with shortest travel time
# - calculate travel time for each arc
# kwargs:
# - levels - the types of highways that will be kept, see ROAD_CLASSES in OSM for more details
# - classSpeeds - dict of the speed (km/hr) of each class of arc
# - classOffRoadAccess - dict indicating which nodes of different arc classes can be used to get on and off road; node will not have off-road access if all connecting arcs have class returning false in dict
# note that distance calculations done here are different (more accurate) from those in JEMSS
function convertOsmNetwork(OsmFilename::String;
	levels::Set{Int} = Set{Int}(), classSpeeds::Dict{Int,Int} = OSM.SPEED_ROADS_RURAL, classOffRoadAccess::Dict{Int,Bool} = Dict{Int,Bool}())
	
	assert(isfile(OsmFilename))
	
	if levels == Set()
		levels = OSM.ROAD_CLASSES |> values |> unique |> sort
	end
	
	# if classOffRoadAccess not given, assume it is true for all levels
	if classOffRoadAccess == Dict()
		classOffRoadAccess = Dict([level => true for level in levels])
	end
	
	# check that input dicts contain all levels
	for level in levels
		assert(haskey(classSpeeds, level))
		assert(haskey(classOffRoadAccess, level))
	end
	
	# convert classSpeeds (km/hr) to classInvSpeeds (days/metre)
	# this is needed as the edge weights (network.w) are in metres
	classInvSpeeds = Dict() # days / metre
	for (class, speed) in classSpeeds
		classInvSpeeds[class] = 1 / (speed * 1000 * 24)
	end
	
	#############################################
	
	## read osm file
	
	(nodesLLA, highways, buildings, features) = OSM.getOSMData(OsmFilename)
	xdoc = OSM.parseMapXML(OsmFilename)
	boundsLLA = OSM.getBounds(xdoc)
	
	# convert from LLA to ENU coordinates
	lla_reference = center(boundsLLA)
	nodesENU = ENU(nodesLLA, lla_reference)
	
	# create graph/network of only specified highway levels
	classes = OSM.roadways(highways)
	network = OSM.createGraph(nodesENU, highways, classes, levels)
	
	graph = network.g # graph type: Graphs.GenericIncidenceList
	assert(graph.is_directed)
	edges = OSM.getEdges(network)
	
	revNetwork = OSM.createGraph(nodesENU, highways, classes, levels, true)
	revGraph = revNetwork.g
	
	#############################################
	
	## find largest strongly connected component in graph, will only keep vertices and edges from this
	
	# create lightGraph from graph
	assert(graph.is_directed)
	numVertices = length(graph.vertices)
	lightGraph = LightGraphs.DiGraph(numVertices)
	# add vertices and edges
	for i = 1:numVertices
		# add edges outgoing from vertex i
		for edge in graph.inclist[i]
			LightGraphs.add_edge!(lightGraph, i, edge.target.index)
		end
	end
	# find connected components; note that Graphs.strongly_connected_components() did not work correctly (the output was not what I expected), so LightGraphs.strongly_connected_components() is needed instead
	components = LightGraphs.strongly_connected_components(lightGraph)
	componentSizes = [length(component) for component in components]
	(largestComponentSize, componentIndex) = findmax(componentSizes)
	largestComponent = components[componentIndex] # largest strongly connected component
	rmComponents = components[vcat(1:componentIndex-1, componentIndex+1:length(components))]
	
	# # show some summary info on strongly connected components
	# @show largestComponentSize
	# @show size(rmComponents) # number of components to remove
	# @show maximum([length(rmComponents[i]) for i = 1:length(rmComponents)]) # size of largest component to remove
	
	# find nodes and arcs to remove
	# remove nodes outside of largestComponent, and remove arcs connected to these nodes
	numVertices = length(graph.vertices)
	numEdges = graph.nedges # length(edges)
	keepVertex = fill(true, numVertices) # keepVertex[i] = true if vertex i is in largestComponent, false otherwise
	keepEdge = fill(true, numEdges) # keepEdge[i] = true if edge i is in largestComponent, false otherwise
	for rmComponent in rmComponents
		for i in rmComponent
			keepVertex[i] = false
			for edge in graph.inclist[i]
				keepEdge[edge.index] = false
			end
			for edge in revGraph.inclist[i]
				keepEdge[edge.index] = false
			end
		end
	end
	
	# remove duplicate arcs in largestComponent, keep arc with shortest travel time
	t(edgeIndex::Int) = network.w[edgeIndex] * classInvSpeeds[network.class[edgeIndex]] # travel time of edge, days
	for i in largestComponent
		vertexEdges = graph.inclist[i]
		n = length(vertexEdges)
		for j = 1:n-1, k = j+1:n
			edge1 = vertexEdges[j]
			if keepEdge[edge1.index]
				edge2 = vertexEdges[k]
				if edge1.target.index == edge2.target.index
					# remove either edge1 or edge2, keep the edge with shorter travel time
					if t(edge1.index) <= t(edge2.index)
						keepEdge[edge2.index] = false
					else
						keepEdge[edge1.index] = false
					end
				end
			end
		end
	end
	
	# keep track of new and old vertex and edge indices
	numKeepVertices = sum(keepVertex) # = length(largestComponent)
	vertexNewIndex = fill(nullIndex, numVertices)
	vertexNewIndex[keepVertex] = 1:numKeepVertices
	# vertexOldIndex = find(keepVertex)
	numKeepEdges = sum(keepEdge)
	edgeNewIndex = fill(nullIndex, numEdges)
	edgeNewIndex[keepEdge] = 1:numKeepEdges
	# edgeOldIndex = find(keepEdge)
	
	#############################################
	
	## convert network into new format, only keeping nodes and edges in largestComponent
	
	numNodes = numKeepVertices
	numArcs = numKeepEdges
	
	nodes = Vector{Node}(numNodes)
	for v in graph.vertices[keepVertex]
		i = vertexNewIndex[v.index]
		assert(i != nullIndex)
		nodes[i] = Node()
		nodes[i].index = i
		nodes[i].location = Location()
		nodes[i].location.x = nodesLLA[v.key].lon
		nodes[i].location.y = nodesLLA[v.key].lat
		nodes[i].offRoadAccess = false # will later change to true where relevant
	end
	
	arcs = Vector{Arc}(numArcs)
	travelTimes = Vector{Float}(numArcs) # travelTimes[i] = travel time along arcs[i]
	for edge in edges[keepEdge]
		i = edgeNewIndex[edge.index]
		assert(i != nullIndex)
		arcs[i] = Arc()
		arcs[i].index = i
		arcs[i].fromNodeIndex = vertexNewIndex[edge.source.index]
		arcs[i].toNodeIndex = vertexNewIndex[edge.target.index]
		assert(arcs[i].fromNodeIndex != nullIndex)
		assert(arcs[i].toNodeIndex != nullIndex)
		
		weight = network.w[edge.index]
		class = network.class[edge.index]
		assert(weight == distance(nodesENU, edge.source.key, edge.target.key))
		travelTimes[i] = weight * classInvSpeeds[class]
		
		if classOffRoadAccess[class]
			nodes[arcs[i].fromNodeIndex].offRoadAccess = true
			nodes[arcs[i].toNodeIndex].offRoadAccess = true
		end
		
		if travelTimes[i] <= 0
			warn("output arc ", i, " has travel time <= 0")
		end
	end
	
	return nodes, arcs, travelTimes # travelTimes has days as units
end

# example use of convertOsmNetwork function
function convertOsmNetworkExample(OsmFilename::String)
	# the 'levels' of highways that will be kept, see ROAD_CLASSES in OSM for more details
	# e.g., levels = Set(1:5) are motorway to tertiary, including links
	# levels = Set(1:6) # roads 1:6 are motorway to residential and unclassified, including links
	levels = Set(1:6)
	
	# speed of each class of highway, in km/hr
	classSpeeds = Dict(
		1 => 105,
		2 => 105,
		3 => 85,
		4 => 65,
		5 => 55,
		6 => 55
	)
	# classSpeeds = OSM.SPEED_ROADS_RURAL # km/hr
	# classSpeeds = OSM.SPEED_ROADS_URBAN # km/hr
	
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
	
	(nodes, arcs, travelTimes) = convertOsmNetwork(OsmFilename;
		levels = levels, classSpeeds = classSpeeds, classOffRoadAccess = classOffRoadAccess)
	
	# # save arcs and nodes files
	# writeArcsFile(arcsFilename, arcs, travelTimes, "directed")
	# writeNodesFile(nodesFilename, nodes)
	
	return nodes, arcs, travelTimes
end
