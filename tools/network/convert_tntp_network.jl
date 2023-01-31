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

# for converting tntp files with nodes and road network to format usable for JEMSS
# tntp files from https://github.com/bstabler/TransportationNetworks

using JEMSS

include("graph_tools.jl")

function readTntpNodesFile(tntpNodesFilename::String; delim::Char='\t')
    # expected node file format (if tab delimited):
    # node	x	y
    # 1	153.2983311	-28.00211673	;
    # 2	153.5293645	-28.19705754	;
    # etc.

    # read nodes file
    println("Reading nodes file")
    data = readDlmFile(tntpNodesFilename; delim=delim)
    cols = findall(.!isempty.(data[1, :])) # need to ignore any empty headers
    table = Table("tntpNodes", data[1, cols], data[2:end, cols])
    numNodes = size(table.data, 1)
    nodes = Vector{Node}(undef, numNodes)
    for i = 1:numNodes
        nodes[i] = Node()
        nodes[i].index = table.columns["node"][i]
        nodes[i].location.x = table.columns["x"][i]
        nodes[i].location.y = table.columns["y"][i]
        @assert(nodes[i].index == i)
    end

    return nodes
end

function readTntpNetworkFile(tntpNetworkFilename::String; delim::Char='\t')
    # expected network file format (if tab delimited):
    # <NUMBER OF ZONES> 1068
    # <NUMBER OF NODES> 4807
    # <FIRST THRU NODE> 1069
    # <NUMBER OF LINKS> 11140
    # [...intermediate lines...]
    # ~	from	to	capacity (per lane)	length (km)	fftt (min)	alpha	beta	ff speed (km/h)	critical speed (km/h)	lanes
    #	 1	1371	900	0.300	0.327	0.282	4	55	42.9	2	;
    #	 2	2012	1600	0.260	0.173	0.282	4	90	70.2	2	;
    # etc.

    # read network file
    println("Reading network file")
    data = readDlmFile(tntpNetworkFilename; delim=delim)

    # first few lines have metadata
    f(r::Regex, s::AbstractString) = Meta.parse(match(r, s).captures[1])
    numZones = f(r"<NUMBER OF ZONES> (\d+)", data[1, 1])
    numNodes = f(r"<NUMBER OF NODES> (\d+)", data[2, 1])
    firstThruNode = f(r"<FIRST THRU NODE> (\d+)", data[3, 1])
    @assert(firstThruNode == numZones + 1)
    numArcs = f(r"<NUMBER OF LINKS> (\d+)", data[4, 1])

    # find and read table data
    row = first(findall(data[:, 1] .== "~"))
    header = data[row, :]
    cols = findall(.!isempty.(header))[2:end] # need to ignore any empty headers, and "~"
    table = Table("tntpNetwork", data[row, cols], data[row+1:end, cols])
    lengthColString = first(filter(s -> occursin("length", s), table.header)) # to match "length", "length (km)", etc
    rowsFields = tableRowsFieldDicts(table, setdiff(table.header, ["from", "to"])) # will keep "length" field
    numArcs = size(table.data, 1)
    arcs = Vector{Arc}(undef, numArcs)
    travelTimes = Vector{Float}(undef, numArcs)
    for i = 1:numArcs
        arcs[i] = Arc()
        arcs[i].index = i
        arcs[i].fromNodeIndex = table.columns["from"][i]
        arcs[i].toNodeIndex = table.columns["to"][i]
        arcs[i].distance = table.columns[lengthColString][i]
        arcs[i].fields = rowsFields[i]
    end

    return arcs, numZones
end

# Process nodes and arcs from tntp, processing performed:
# - remove all zones and their connected arcs
# - keep the largest strongly connected component
function convertTntpNetwork!(nodes::Vector{Node}, arcs::Vector{Arc}, numZones::Int)

    @assert(!graphContainsDuplicateArcs(arcs))

    # remove all zones and their connected arcs
    nodeFilter(node::Node) = node.index > numZones
    arcFilter(arc::Arc) = arc.fromNodeIndex > numZones && arc.toNodeIndex > numZones
    graphRemoveElts!(nodes, arcs; nodeFilter=nodeFilter, arcFilter=arcFilter)

    graphRemoveDisconnectedArcs!(nodes, arcs) # just in case, but should not be necessary
    graphKeepLargestComponent!(nodes, arcs)
end

function convertTntpNetworkFile(tntpNodesFilename::String, tntpNetworkFilename::String;
    delim::Char='\t')
    nodes = readTntpNodesFile(tntpNodesFilename; delim=delim)
    (arcs, numZones) = readTntpNetworkFile(tntpNetworkFilename; delim=delim)
    convertTntpNetwork!(nodes, arcs, numZones)
    return nodes, arcs
end
