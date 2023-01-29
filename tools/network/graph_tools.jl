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

# tools to help manipulate graphs

using JEMSS
import LightGraphs
using SparseArrays

travelModeString(i::Int) = "mode_$i"

function checkNodeIndices(nodes::Vector{Node})
    result = all(i -> i == nodes[i].index, 1:length(nodes)) # true if ok
    result || @warn("Node indices are not 1:n")
    return result
end

function checkArcIndices(arcs::Vector{Arc})
    result = all(i -> i == arcs[i].index, 1:length(arcs)) # true if ok
    result || @warn("Arc indices are not 1:n")
    return result
end

# check arc.fromNodeIndex and arc.toNodeIndex values are in expected range
function checkArcNodeIndices(numNodes::Int, arcs::Vector{Arc})
    result = all(i -> (1 <= arcs[i].fromNodeIndex <= numNodes) && (1 <= arcs[i].toNodeIndex <= numNodes), 1:length(arcs)) # true if ok
    result || @warn("Arc from/to node indices are not all within 1:numNodes")
    return result
end
checkArcNodeIndices(nodes::Vector{Node}, arcs::Vector{Arc}) = checkArcNodeIndices(length(nodes), arcs)

# get number of travel modes for an arc
function getNumTravelModes(arc::Arc)
    i = 1
    while haskey(arc.fields, travelModeString(i))
        i += 1
    end
    numModes = i - 1
    @assert(numModes >= 1, "Arcs should contain at least one travel mode field, e.g. '$(travelModeString(1))'")
    return numModes
end
getNumTravelModes(arcs::Vector{Arc}) = getNumTravelModes(arcs[1]) # assume all arcs have same fields

# get travel times for each travel mode of an arc
function getArcTravelTimes(arc::Arc)
    numModes = getNumTravelModes(arcs)
    travelTimes = Vector{Float}(undef, numModes)
    for i = 1:numModes
        travelTimes[i] = arc.fields[travelModeString(i)]
    end
    return travelTimes
end

# return vector of travel times for arcs,
# first dimension is for travel modes, second is for arcs
function getArcsTravelTimes(arcs::Vector{Arc})
    numArcs = length(arcs)
    numModes = getNumTravelModes(arcs)
    travelTimes = Array{Float,2}(undef, numModes, numArcs)
    for i = 1:numModes
        modeString = travelModeString(i)
        for j = 1:numArcs
            travelTimes[i, j] = arcs[j].fields[modeString]
        end
    end
    return travelTimes
end

function nodesAddField!(nodes::Vector{Node}, field::String; default::Any=nothing, overwrite::Bool=false)
    if !overwrite && haskey(nodes[1].fields, field)
        return
    end
    for node in nodes
        node.fields[field] = default
    end
end

function arcsAddField!(arcs::Vector{Arc}, field::String; default::Any=nothing, overwrite::Bool=false)
    if !overwrite && haskey(arcs[1].fields, field)
        return
    end
    for arc in arcs
        arc.fields[field] = default
    end
end

nodesDeleteField!(nodes::Vector{Node}, field::String) =
    for node in nodes
        delete!(node.fields, field)
    end
arcsDeleteField!(arcs::Vector{Arc}, field::String) =
    for arc in arcs
        delete!(arc.fields, field)
    end

function JEMSS.writeArcsFile(arcsFilename::String, arcs::Vector{Arc}, arcForm::String)
    travelTimes = getArcsTravelTimes(arcs)
    writeArcsFile(arcsOutputFilename, arcs, travelTimes, arcForm)
end

function createLightGraph(numNodes::Int, arcs::Vector{Arc})
    lightGraph = LightGraphs.DiGraph(numNodes)
    for arc in arcs
        LightGraphs.add_edge!(lightGraph, arc.fromNodeIndex, arc.toNodeIndex)
    end
    return lightGraph
end
createLightGraph(nodes::Vector{Node}, arcs::Vector{Arc}) = createLightGraph(length(nodes), arcs)

# tag the elements (nodes and arcs) that contribute to the largest strongly connected component
# tagging is done by adding a field with name given by kwarg `lsccHeader`
function graphTagLargestComponentElts!(nodes::Vector{Node}, arcs::Vector{Arc};
    lsccHeader::String="in_largest_component")

    # shorthand
    numNodes = length(nodes)
    numArcs = length(arcs)

    lightGraph = createLightGraph(nodes, arcs)

    # find strongly connected components; note that Graphs.strongly_connected_components() did not work correctly (the output was not what I expected), so LightGraphs.strongly_connected_components() is needed instead
    components = LightGraphs.strongly_connected_components(lightGraph)
    largestComponent = components[argmax([length(c) for c in components])] # node indices of largest strongly connected component

    # tag nodes and arcs in largest strongly connected component
    tagNode = fill(false, numNodes) # tagNode[i] = true if nodes[i] in largest component
    tagArc = fill(false, numArcs) # tagArc[i] = true if arcs[i] in largest component
    for i in largestComponent
        tagNode[i] = true
    end
    for (i, arc) in enumerate(arcs)
        if tagNode[arc.fromNodeIndex] && tagNode[arc.toNodeIndex]
            tagArc[i] = true
        end
    end

    # tag by adding fields to nodes and arcs
    for i = 1:numNodes
        nodes[i].fields[lsccHeader] = tagNode[i]
    end
    for i = 1:numArcs
        arcs[i].fields[lsccHeader] = tagArc[i]
    end

    return
end

# keep nodes and arcs where filters return true
# indices will be renumbered
# arcs with either (to/from) node removed are not automatically removed
function graphRemoveElts!(nodes::Vector{Node}, arcs::Vector{Arc};
    nodeFilter::Function=(x -> true), arcFilter::Function=(x -> true))

    # before filtering, create mapping from old indices to new
    keepNodeIndices = findall(nodeFilter, nodes)
    nodeNewIndex = fill(nullIndex, length(nodes))
    nodeNewIndex[keepNodeIndices] = 1:length(keepNodeIndices) # nodeNewIndex[i] gives new node index for old nodes[i]

    # filter nodes, renumber nodes from 1 to n
    filter!(nodeFilter, nodes)
    for (i, node) in enumerate(nodes)
        node.index = i
    end

    # filter arcs, renumber nodes from 1 to n
    filter!(arcFilter, arcs)
    for (i, arc) in enumerate(arcs)
        arc.index = i
    end

    # connect arcs to nodes
    # arc.fromNodeIndex or arc.toNodeIndex will become nullIndex if node does not exist
    for arc in arcs
        arc.fromNodeIndex = nodeNewIndex[arc.fromNodeIndex]
        arc.toNodeIndex = nodeNewIndex[arc.toNodeIndex]
    end

    # check that all needed nodes exist for filtered arcs
    checkArcNodeIndices(nodes, arcs)

    return
end

# remove elements outside of largest strongly connected component
function graphKeepLargestComponent!(nodes::Vector{Node}, arcs::Vector{Arc};
    lsccHeader::String="in_largest_component")
    graphTagLargestComponentElts!(nodes, arcs; lsccHeader=lsccHeader)
    nodeFilter(node::Node) = node.fields[lsccHeader]
    arcFilter(arc::Arc) = arc.fields[lsccHeader]
    graphRemoveElts!(nodes, arcs, nodeFilter=nodeFilter, arcFilter=arcFilter)
    nodesDeleteField!(nodes, lsccHeader)
    arcsDeleteField!(arcs, lsccHeader)
    return
end

# remove arcs from graph that are connected to/from nodes outside of range
function graphRemoveDisconnectedArcs!(nodes::Vector{Node}, arcs::Vector{Arc})
    numNodes = length(nodes)
    arcFilter(arc::Arc) = (1 <= arc.fromNodeIndex <= numNodes) && (1 <= arc.toNodeIndex <= numNodes)
    graphRemoveElts!(nodes, arcs, arcFilter=arcFilter)
end

# return dict with tuple of node indices (from, to) as keys, and vector of arc indices as values
function graphNodePairArcIndices(arcs::Vector{Arc})
    nodePairArcIndices = Dict{Tuple{Int,Int},Vector{Int}}() # nodePairArcIndices[i,j] gives indices of arcs going from node i to node j
    for arc in arcs
        i = arc.fromNodeIndex
        j = arc.toNodeIndex
        if !haskey(nodePairArcIndices, (i, j))
            nodePairArcIndices[(i, j)] = []
        end
        push!(nodePairArcIndices[(i, j)], arc.index)
    end
    return nodePairArcIndices
end

# find duplicate arcs (same 'to' and 'from' nodes),
# return vector of vector of indices of duplicate arcs
function graphFindDuplicateArcs(arcs::Vector{Arc})
    nodePairArcIndices = graphNodePairArcIndices(arcs)
    duplicateArcsSets = []
    for (nodePair, arcIndices) in nodePairArcIndices
        if length(arcIndices) > 1
            push!(duplicateArcsSets, arcIndices)
        end
    end
    return duplicateArcsSets
end

graphContainsDuplicateArcs(arcs::Vector{Arc}) = graphFindDuplicateArcs(arcs) != []
function checkGraphContainsNoDuplicateArcs(arcs::Vector{Arc})
    result = !graphContainsDuplicateArcs(arcs) # true if ok
    result || @warn("There are duplicate arcs.")
    return result
end

# merge arc2 into arc1, keeping minimum travel time (from both arcs) for each travel mode in arc1
# arc2 will need to be removed in separate step
function mergeDuplicateArcs!(arc1::Arc, arc2::Arc)
    @assert(arc1.fromNodeIndex == arc2.fromNodeIndex)
    @assert(arc1.toNodeIndex == arc2.toNodeIndex)
    numModes = getNumTravelModes(arc1)
    for i = 1:numModes
        modeString = travelModeString(i)
        arc1.fields[modeString] = min(arc1.fields[modeString], arc2.fields[modeString])
    end
    # assume that arc.distance is same for both arcs
end
mergeDuplicateArcs!(arcs::Vector{Arc}, i::Int, j::Int) = mergeDuplicateArcs!(arcs[i], arcs[j])

# merge any duplicate arcs into one
# will take the shortest travel time from any merged arcs, for each travel mode
function graphMergeDuplicateArcs!(nodes::Vector{Node}, arcs::Vector{Arc};
    mergeArcsFunction::Function=mergeDuplicateArcs!, mergeHeader::String="merge_result")

    arcsAddField!(arcs, mergeHeader, default="unmerged")

    # find which arcs are duplicates
    duplicateArcsSets = graphFindDuplicateArcs(arcs) # duplicateArcsSets[i] has a vector of indices of duplicate arcs
    for duplicateArcs in duplicateArcsSets
        i = duplicateArcs[1]
        arcs[i].fields[mergeHeader] = "merged"
        for j in duplicateArcs[2:end]
            mergeArcsFunction(arcs, i, j)
            arcs[j].fields[mergeHeader] = "merged_in"
        end
    end

    # remove duplicate arcs
    arcFilter(arc::Arc) = arc.fields[mergeHeader] != "merged_in"
    graphRemoveElts!(nodes, arcs, arcFilter=arcFilter)
end

# Divide an arc into even sections a given number of times, 0 divides makes no change.
# Only some nodes and arcs fields are updated (node offRoadAccess, arc travel times and osm weight), the rest are unchanged but may need changing.
function graphDivideArc!(nodes::Vector{Node}, arcs::Vector{Arc}; arcIndex::Int=0, numDivides::Int=0)
    if numDivides <= 0
        return
    end
    @assert(1 <= arcIndex <= length(arcs))

    # shorthand
    numNodes = length(nodes)
    numArcs = length(arcs)
    arc = arcs[arcIndex]
    fromNode = nodes[arc.fromNodeIndex]
    toNode = nodes[arc.toNodeIndex]

    # keep track of which nodes have been added, and which arcs have been divided
    addedHeader = "added"
    divArcHeader = "divided_arc_index"
    nodesAddField!(nodes, addedHeader; default=false)
    arcsAddField!(arcs, divArcHeader; default=nullIndex)

    # change some fields of arc
    for i = 1:getNumTravelModes(arc)
        arc.fields[travelModeString(i)] /= numDivides + 1
    end
    if haskey(arc.fields, "osm_weight")
        arc.fields["osm_weight"] /= numDivides + 1
    end
    arc.distance /= numDivides + 1
    arc.fields[divArcHeader] = arc.index

    tempNode = deepcopy(fromNode) # could copy any node
    tempNode.fields[addedHeader] = true
    tempArc = deepcopy(arc)
    # some fields of these new nodes and arcs will have incorrect/irrelevant information and should be ignored

    for i = 1:numDivides
        # add node
        newNode = deepcopy(tempNode)
        newNode.index = numNodes + i
        newNode.location = linearInterpLocation(fromNode.location, toNode.location, 0.0, 1.0, i / (numDivides + 1))
        newNode.offRoadAccess = fromNode.offRoadAccess && toNode.offRoadAccess
        push!(nodes, newNode)

        # add arc
        newArc = deepcopy(tempArc)
        newArc.index = numArcs + i
        newArc.fromNodeIndex = newNode.index
        newArc.toNodeIndex = newNode.index + 1
        push!(arcs, newArc)
    end
    arc.toNodeIndex = numNodes + 1
    arcs[end].toNodeIndex = toNode.index
end

# Divide arcs into multiple arcs (with added nodes) so that the the maximum arc travel time <= maxArcTravelTime.
function graphDivideArcs!(nodes::Vector{Node}, arcs::Vector{Arc}; maxArcTravelTime::Float=Inf)
    if maxArcTravelTime == Inf
        return
    end
    @assert(maxArcTravelTime > 0)

    nodePairArcIndex = Dict{Tuple{Int,Int},Int}()
    for arc in arcs
        @assert(!haskey(nodePairArcIndex, (arc.fromNodeIndex, arc.toNodeIndex)), arc) # duplicate arcs not allowed
        nodePairArcIndex[arc.fromNodeIndex, arc.toNodeIndex] = arc.index
    end

    travelTimes = getArcsTravelTimes(arcs) # values before dividing arcs, actual values will change (but not in this var)
    arcsDivided = fill(false, length(arcs))
    for (node1, node2) in keys(nodePairArcIndex)
        if arcsDivided[nodePairArcIndex[node1, node2]]
            continue
        end # arc has already been divided
        numNodes = length(nodes) # keep track of number of nodes before current divide
        if !haskey(nodePairArcIndex, (node2, node1)) # no arc in opposite direction
            i = nodePairArcIndex[node1, node2]
            numDivides = ceil(Int, maximum(travelTimes[:, i]) / maxArcTravelTime) - 1
            if numDivides < 1
                continue
            end
            graphDivideArc!(nodes, arcs; arcIndex=i, numDivides=numDivides)
            arcsDivided[i] = true
        else
            # there are two arcs in opposite directions, need to split both simultaneously
            arc1 = nodePairArcIndex[node1, node2]
            arc2 = nodePairArcIndex[node2, node1]
            numDivides = ceil(Int, maximum(travelTimes[:, [arc1, arc2]]) / maxArcTravelTime) - 1 # look at maximum travel time of both arcs
            if numDivides < 1
                continue
            end
            graphDivideArc!(nodes, arcs; arcIndex=arc1, numDivides=numDivides)
            graphDivideArc!(nodes, arcs; arcIndex=arc2, numDivides=numDivides)
            arcsDivided[arc1] = true
            arcsDivided[arc2] = true

            # have created extra nodes when dividing arc2; need to remove these and instead use nodes that were created when dividing arc1
            for i = 1:numDivides
                pop!(nodes)
            end
            for i = 1:numDivides
                arc = arcs[end-i+1]
                arc.fromNodeIndex = numNodes + i
                i == 1 || (arc.toNodeIndex = numNodes + i - 1)
            end
            arcs[arc2].toNodeIndex = numNodes + numDivides
        end
    end

    # check max arc travel time
    travelTimes = getArcsTravelTimes(arcs)
    @assert(maximum(travelTimes) <= maxArcTravelTime)

    checkNodeIndices(nodes)
    checkArcIndices(arcs)
    checkArcNodeIndices(nodes, arcs)
    checkGraphContainsNoDuplicateArcs(arcs)
end

# Tag all nodes and arcs that are used in any shortest path:
# from an origin node to a chosen node,
# from a chosen node to a destination node, or
# from an origin node to a destination node
# This assumes that any duplicate arcs have been removed.
# This version requires LightGraphs.dijkstra_shortest_paths function to allow distmx::Dict, for sake of reducing run time
function graphTagSpElts!(nodes::Vector{Node}, arcs::Vector{Arc};
    spHeader::String="in_a_sp", chosenNodes::Vector{Int}=Int[], originNodes::Vector{Int}=Int[], destNodes::Vector{Int}=Int[], dijkstraSpDistmxAllowDict::Bool=false)

    @assert(!graphContainsDuplicateArcs(arcs), "Cannot handle duplicate arcs")
    if dijkstraSpDistmxAllowDict
        @info("Setting `dijkstraSpDistmxAllowDict` = true requires LightGraphs.dijkstra_shortest_paths function to be changed, to allow arg `distmx` to be type Dict")
    end

    # shorthand
    numArcs = length(arcs)
    numNodes = length(nodes)

    # keep track of all tagged nodes and arcs
    nodeTagged = fill(false, numNodes)
    arcTagged = fill(false, numArcs)

    # to be used, and reset, for each SP tree:
    tagNode = Vector{Bool}(undef, numNodes)
    tagArc = Vector{Bool}(undef, numArcs)

    nodePairArcIndex = Dict{Tuple{Int,Int},Int}()
    revNodePairArcIndex = Dict{Tuple{Int,Int},Int}()
    for arc in arcs
        nodePairArcIndex[arc.fromNodeIndex, arc.toNodeIndex] = arc.index
        revNodePairArcIndex[arc.toNodeIndex, arc.fromNodeIndex] = arc.index
    end

    # # using sparse matrices takes much longer to populate than using dict
    # nodePairArcIndex = spzeros(Int, numNodes, numNodes)
    # for arc in arcs
    # nodePairArcIndex[arc.fromNodeIndex, arc.toNodeIndex] = arc.index
    # end
    # revNodePairArcIndex = transpose(nodePairArcIndex)

    lightGraph = createLightGraph(nodes, arcs)
    revLightGraph = LightGraphs.reverse(lightGraph)

    travelTimes = getArcsTravelTimes(arcs)
    numModes = size(travelTimes, 1)

    # create distance matrix, to be reused
    println("Setting arc time matrix / dict")
    if dijkstraSpDistmxAllowDict
        # using this with LightGraphs requires dijkstra_shortest_paths function distmx arg type to be changed
        arcTimes = Dict{Tuple{Int,Int},Float}()
        revArcTimes = Dict{Tuple{Int,Int},Float}()
        for arc in arcs
            if mod(arc.index, 1000) == 0
                print("\r\tarc: ", arc.index)
            end
            arcTimes[arc.fromNodeIndex, arc.toNodeIndex] = nullTime # to be set later
            revArcTimes[arc.toNodeIndex, arc.fromNodeIndex] = nullTime # to be set later
        end
        println()
    else
        # use sparse matrix
        # this is much slower to populate than a dict, but works with unedited version of LightGraphs.dijkstra_shortest_paths
        arcTimes = spzeros(Float, numNodes, numNodes)
        revArcTimes = spzeros(Float, numNodes, numNodes)
        for arc in arcs
            if mod(arc.index, 1000) == 0
                print("\r\tarc: ", arc.index)
            end
            arcTimes[arc.fromNodeIndex, arc.toNodeIndex] = nullTime # to be set later
            revArcTimes[arc.toNodeIndex, arc.fromNodeIndex] = nullTime # to be set later
        end
        println()
    end

    # Tag all nodes and arcs (by setting vars tagNode and tagArc) that are in a shortest path from the originNode to any node in destNodeSets
    # To work with destination node (instead of origin), provide a reversed lightGraph and transposed arcTimes and nodePairArcIndex
    function graphTagSpTreeElts!(lightGraph::LightGraphs.DiGraph{Int}, originNode::Int, destNodeSets::Vector{Vector{Int}}, arcTimes::Union{SparseMatrixCSC{Float,Int},Dict{Tuple{Int,Int},Float}}, nodePairArcIndex::Dict{Tuple{Int,Int},Int})
        # global tagNode, tagArc
        i = originNode # shorthand
        spData = LightGraphs.dijkstra_shortest_paths(lightGraph, i, arcTimes)
        spPreds = spData.parents # spPreds[j] gives index of node before j on shortest path from i to j
        tagNode[:] .= false
        tagNode[i] = true
        tagArc[:] .= false
        for nodeIndices in destNodeSets
            for j in nodeIndices
                while !tagNode[j]
                    tagNode[j] = true
                    k = j
                    j = spPreds[j]
                    tagArc[nodePairArcIndex[j, k]] = true
                end
            end
        end
        nodeTagged .|= tagNode
        arcTagged .|= tagArc
        return
    end

    println("Calculating shortest paths")
    for modeIndex in 1:numModes

        println("\tmode: $modeIndex")

        # set distance matrix values
        for arc in arcs
            arcTimes[arc.fromNodeIndex, arc.toNodeIndex] = travelTimes[modeIndex, arc.index]
            revArcTimes[arc.toNodeIndex, arc.fromNodeIndex] = travelTimes[modeIndex, arc.index]
        end

        for (i, originNode) in enumerate(originNodes)
            print("\r\torigin node: $i (of $(length(originNodes)))")
            # tag nodes in a shortest path from node originNode to any chosen/dest nodes
            graphTagSpTreeElts!(lightGraph, originNode, [chosenNodes, destNodes], arcTimes, nodePairArcIndex)
        end
        println()

        for (i, destNode) in enumerate(destNodes)
            print("\r\tdest node: $i (of $(length(destNodes)))")
            # tag nodes in a shortest path from any chosen/origin node to node destNode
            graphTagSpTreeElts!(revLightGraph, destNode, [chosenNodes, originNodes], revArcTimes, revNodePairArcIndex)
        end
        println()

    end

    # use vars `nodeTagged` and `arcTagged` to set value of new field for nodes and arcs
    for i = 1:numNodes
        nodes[i].fields[spHeader] = nodeTagged[i]
    end
    for i = 1:numArcs
        arcs[i].fields[spHeader] = arcTagged[i]
    end

    return
end
