# Create nodes and arcs files.

using JEMSS
using Shapefile
using GeoInterface

dir = @__DIR__

jemssDir = JEMSS.jemssDir # shorthand
include(joinpath(jemssDir, "tools/network/convert_osm_network.jl"))

## parameters

osmFilename = joinpath(jemssDir, "data/cities/edmonton/data/roads/osm 2019/edmonton_main_roads_tertiary+.osm")

levels = Set(1:6) # osm classes 1:8 : motorway, trunk, primary, secondary, tertiary, unclassified/residential, service, pedestrian/living street

# arc speed (km/hr) for each travel mode and each arc class
v2 = [102, 68, 49, 45, 45, 45] # osm classes 1:8 : motorway, trunk, primary, secondary, tertiary, unclassified/residential, service, pedestrian/living street
v2 *= 0.7 # reduce speeds
v = [
    round.(Int, v2 / 0.7), # travel mode 1
    round.(Int, v2) # travel mode 2
]
classSpeeds = [Dict([j => Float(v[i][j]) for j = 1:length(v[i])]) for i = 1:length(v)]

# tag which nodes can be used for off-road access
v = [false, false, true, true, true, true, true, true] # v[i] is true if node at either end of arc with class = i can be used for off-road access, false otherwise
classOffRoadAccess = Dict([i => v[i] for i = 1:length(v)])

# only keep nodes and arcs that provide a shortest path between a hospital/station and arcs of given classes
spTagClass = [true, true, true, true, true, false, false, false]

hospitalsFilename = joinpath(jemssDir, "data/cities/edmonton/models/1/hospitals/hospitals_1.csv")
stationsFilename = joinpath(jemssDir, "data/cities/edmonton/models/1/stations/stations_1.csv")
mapFilename = joinpath(jemssDir, "data/cities/edmonton/models/1/maps/map_1.csv")

# shapefile of region border
borderFilename = joinpath(jemssDir, "data/cities/edmonton/data/border/census 2011/edmonton.shp")

# divide arcs that are long into multiple arcs, according to maxArcTravelTime value
maxArcTravelTime = 10 / 60 / 60 / 24 # convert seconds to days

## processing

(nodes, arcs) = readOsmNetworkFile(osmFilename; levels=levels);
convertOsmNetwork!(nodes, arcs; levels=levels, classSpeeds=classSpeeds, classOffRoadAccess=classOffRoadAccess)

print("Dividing arcs")
graphDivideArcs!(nodes, arcs; maxArcTravelTime=maxArcTravelTime)
println(" ... done.")

# get nodes that are nearest to hospitals and stations
begin # get od nodes
    hospitals = readHospitalsFile(hospitalsFilename)
    stations = readStationsFile(stationsFilename)
    mp = readMapFile(mapFilename)
    odNodes = Int[]
    for items in [hospitals, stations]
        for item in items
            i, dist = findNearestNode(mp, nodes, item.location) # if this is slow, instead use version of findNearestNode that uses grid
            push!(odNodes, i)
        end
    end
    odNodes = odNodes |> unique |> sort
end # get od nodes

nodeChosen = fill(false, length(nodes));
for arc in arcs
    if spTagClass[arc.fields["osm_class"]]
        nodeChosen[arc.fromNodeIndex] = true
        nodeChosen[arc.toNodeIndex] = true
    end
end

# remove nodes outside region from chosen nodes
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
    if nodeChosen[i]
        location = nodes[i].location
        nodeChosen[i] = pointInRings([location.x, location.y], rings)
    end
end
println("\tdone")
@warn("Code to check that nodes are inside the border is not guaranteed to be correct; please check that output is as expected.")

for odNode in odNodes
    nodeChosen[odNode] = true
end
chosenNodes = findall(nodeChosen);

graphTagSpElts!(nodes, arcs, chosenNodes=chosenNodes, originNodes=odNodes, destNodes=odNodes, dijkstraSpDistmxAllowDict=true)

# only keep elements used in shortest paths
# and arcs with required osm class
nodeFilter(node::Node) = node.fields["in_a_sp"]
arcFilter(arc::Arc) = spTagClass[arc.fields["osm_class"]] || arc.fields["in_a_sp"]
graphRemoveElts!(nodes, arcs, nodeFilter=nodeFilter, arcFilter=arcFilter)

# remove arcs not connected to valid nodes,
# needed as we only kept nodes inside the given region (see border shapefile)
graphRemoveDisconnectedArcs!(nodes, arcs)

# change osm_key field of node to string, to preserve value when writing to file
for node in nodes
    node.fields["osm_key"] = string(node.fields["osm_key"])
end

# remove some fields before saving, to reduce file size and reduce time when reading file
travelTimes = getArcsTravelTimes(arcs) # need to extract travel times before deleting them from arcs
for field in keys(nodes[1].fields)
    nodesDeleteField!(nodes, field)
end
for field in keys(arcs[1].fields)
    arcsDeleteField!(arcs, field)
end

# save nodes and arcs to file
nodesOutputFilename = joinpath(dir, "nodes.csv")
arcsOutputFilename = joinpath(dir, "arcs.csv")
writeNodesFile(nodesOutputFilename, nodes)
writeArcsFile(arcsOutputFilename, arcs, travelTimes, "directed")
