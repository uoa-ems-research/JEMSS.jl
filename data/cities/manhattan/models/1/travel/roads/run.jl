# Create nodes and arcs files.

using JEMSS
using Shapefile
using GeoInterface

dir = @__DIR__

jemssDir = JEMSS.jemssDir # shorthand
include(joinpath(jemssDir, "tools/network/convert_osm_network.jl"))

## parameters

osmFilename = joinpath(jemssDir, "data/cities/manhattan/data/roads/osm 2019/manhattan_main_roads_tertiary+_unclass_resid.osm")

levels = Set(1:6) # osm classes 1:8 : motorway, trunk, primary, secondary, tertiary, unclassified/residential, service, pedestrian/living street

# arc speed (km/hr) for each travel mode and each arc class
v2 = [72, 35, 25, 22, 11, 11] # travel mode 2, classes 1:6
v = [v2 * 4 / 3, v2] # travel modes 1 and 2
classSpeeds = [Dict([j => Float(v[i][j]) for j = 1:length(v[i])]) for i = 1:length(v)]

# tag which nodes can be used for off-road access
v = [false, false, true, true, true, true, true, true] # v[i] is true if node at either end of arc with class = i can be used for off-road access, false otherwise
classOffRoadAccess = Dict([i => v[i] for i = 1:length(v)])

# shapefile of region border
borderFilename = joinpath(jemssDir, "data/cities/manhattan/data/border/census 2010/manhattan_border.shp")

## processing

(nodes, arcs) = readOsmNetworkFile(osmFilename; levels=levels);
convertOsmNetwork!(nodes, arcs; levels=levels, classSpeeds=classSpeeds, classOffRoadAccess=classOffRoadAccess)

# remove nodes outside region from chosen nodes
nodeInBorder = fill(false, length(nodes));
shp = open(borderFilename) do data
    read(data, Shapefile.Handle)
end
@assert(length(shp.shapes) == 1)
coords = GeoInterface.coordinates(shp.shapes[1]) # idk what to call this var
rings = unique(vcat(coords...)) # there can be some repeated rings (not sure why), these need to be removed
# check if a point is in given rings (should be in an odd number of rings, to be on land)
function pointInRings(pt::Vector{Float}, rings::Vector{Vector{Vector{Float}}})::Bool
    count = 0
    for ring in rings
        count += Shapefile.inring(pt, ring)
    end
    return isodd(count)
end
println("Checking which nodes (of $(length(nodes))) are inside border")
for i = 1:length(nodes)
    mod(i, 1000) == 0 ? print("\r node: $i") : nothing
    location = nodes[i].location
    nodeInBorder[i] = pointInRings([location.x, location.y], rings)
end
println("\tdone")
@warn("Code to check that nodes are inside the border is not guaranteed to work; please check that output is as expected.")

nodeFilter(node::Node) = nodeInBorder[node.index]
graphRemoveElts!(nodes, arcs, nodeFilter=nodeFilter)

# remove arcs not connected to valid nodes,
# needed as we only kept nodes inside the given region (see border shapefile)
graphRemoveDisconnectedArcs!(nodes, arcs)

graphKeepLargestComponent!(nodes, arcs)

# # change osm_key field of node to string, to preserve value when writing to file
# for node in nodes node.fields["osm_key"] = string(node.fields["osm_key"]) end

# remove some fields before saving
for field in ["osm_key"]
    nodesDeleteField!(nodes, field)
end
for field in ["osm_weight", "osm_class", "merge_result"]
    arcsDeleteField!(arcs, field)
end

# save nodes and arcs to file
nodesOutputFilename = joinpath(dir, "nodes.csv")
arcsOutputFilename = joinpath(dir, "arcs.csv")
writeNodesFile(nodesOutputFilename, nodes)
writeArcsFile(arcsOutputFilename, arcs, "directed")
