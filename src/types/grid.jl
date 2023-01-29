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

# returns grid indices for given location
function locationToGridIndex(map::Map, grid::Grid, location::Location)
    # fails if location.x == map.xMin or location.y == map.yMin, can fix this later if needed
    ix = Int(ceil((location.x - map.xMin) / map.xRange * grid.nx))
    iy = Int(ceil((location.y - map.yMin) / map.yRange * grid.ny))
    if location.x == map.xMin
        ix = 1
    end
    if location.y == map.yMin
        iy = 1
    end
    @assert(1 <= ix && ix <= grid.nx && 1 <= iy && iy <= grid.ny)
    return ix, iy
end

# place nodes in rectangles of grid
# can specify whether nodes require offRoadAccess = true, in order to be added
function gridPlaceNodes!(map::Map, grid::Grid, nodes::Vector{Node};
    offRoadAccessRequired::Bool=true)
    for node in nodes
        if node.offRoadAccess || !offRoadAccessRequired
            (ix, iy) = locationToGridIndex(map, grid, node.location)
            push!(grid.rects[ix, iy].nodeIndices, node.index)
        end
    end
end

# searches grid for nearest node to given location,
# starting in grid rectangle of location, then
# extending search outward in a rectangle, until
# borders of search rectangle are further from given
# location than the nearest node found
function findNearestNode(map::Map, grid::Grid, nodes::Vector{Node}, location::Location)

    # indices of grid rectangle for given location
    (ix, iy) = locationToGridIndex(map, grid, location)

    # initialise grid search rectangle
    gsr = grid.searchRect
    gsr.ixSearched[1] = ix
    gsr.ixSearched[2] = ix
    gsr.iySearched[1] = iy
    gsr.iySearched[2] = iy
    gsr.ixSearch[1] = ix
    gsr.ixSearch[2] = ix
    gsr.iySearch[1] = iy
    gsr.iySearch[2] = iy
    gsr.xDist[1] = (location.x - (map.xMin + grid.xRange * (ix - 1))) * map.xScale
    gsr.xDist[2] = grid.xRangeDist - gsr.xDist[1]
    gsr.yDist[1] = (location.y - (map.yMin + grid.yRange * (iy - 1))) * map.yScale
    gsr.yDist[2] = grid.yRangeDist - gsr.yDist[1]

    nearestNode = nullIndex
    bestDist = Inf
    node = nullIndex # init
    dist = nullDist # init

    found = false
    skipSearch = false
    while !found

        if !skipSearch
            for i = gsr.ixSearch[1]:gsr.ixSearch[2]
                for j = gsr.iySearch[1]:gsr.iySearch[2]
                    # find nearest node in grid rect
                    nearestNodeInGrid = nullIndex
                    bestSqrDist = Inf
                    sqrDist = 0.0 # init
                    for nodeIndex in grid.rects[i, j].nodeIndices
                        sqrDist = squareDist(map, nodes[nodeIndex].location, location)
                        if sqrDist < bestSqrDist
                            bestSqrDist = sqrDist
                            nearestNodeInGrid = nodeIndex
                        end
                    end
                    dist = sqrt(bestSqrDist)

                    if dist < bestDist
                        bestDist = dist
                        nearestNode = nearestNodeInGrid
                    end
                end
            end
        end
        nearestBorderDist = min(gsr.xDist[1], gsr.xDist[2], gsr.yDist[1], gsr.yDist[2])
        if bestDist <= nearestBorderDist
            found = true
        else
            # extend rectangle search border
            skipSearch = false
            if gsr.xDist[1] == nearestBorderDist
                if gsr.ixSearched[1] > 1 # extend border
                    gsr.ixSearched[1] -= 1
                    gsr.xDist[1] += grid.xRangeDist
                    gsr.ixSearch[1] = gsr.ixSearched[1]
                    gsr.ixSearch[2] = gsr.ixSearched[1]
                    gsr.iySearch[1] = gsr.iySearched[1]
                    gsr.iySearch[2] = gsr.iySearched[2]
                else # cannot extend border
                    gsr.xDist[1] = Inf
                    skipSearch = true
                end
            elseif gsr.xDist[2] == nearestBorderDist
                if gsr.ixSearched[2] < grid.nx # extend border
                    gsr.ixSearched[2] += 1
                    gsr.xDist[2] += grid.xRangeDist
                    gsr.ixSearch[1] = gsr.ixSearched[2]
                    gsr.ixSearch[2] = gsr.ixSearched[2]
                    gsr.iySearch[1] = gsr.iySearched[1]
                    gsr.iySearch[2] = gsr.iySearched[2]
                else # cannot extend border
                    gsr.xDist[2] = Inf
                    skipSearch = true
                end
            elseif gsr.yDist[1] == nearestBorderDist
                if gsr.iySearched[1] > 1 # extend border
                    gsr.iySearched[1] -= 1
                    gsr.yDist[1] += grid.yRangeDist
                    gsr.ixSearch[1] = gsr.ixSearched[1]
                    gsr.ixSearch[2] = gsr.ixSearched[2]
                    gsr.iySearch[1] = gsr.iySearched[1]
                    gsr.iySearch[2] = gsr.iySearched[1]
                else # cannot extend border
                    gsr.yDist[1] = Inf
                    skipSearch = true
                end
            elseif gsr.yDist[2] == nearestBorderDist
                if gsr.iySearched[2] < grid.ny # extend border
                    gsr.iySearched[2] += 1
                    gsr.yDist[2] += grid.yRangeDist
                    gsr.ixSearch[1] = gsr.ixSearched[1]
                    gsr.ixSearch[2] = gsr.ixSearched[2]
                    gsr.iySearch[1] = gsr.iySearched[2]
                    gsr.iySearch[2] = gsr.iySearched[2]
                else # cannot extend border
                    gsr.yDist[2] = Inf
                    skipSearch = true
                end
            end
        end

    end

    return nearestNode, bestDist # (gsr.ixSearched[2] - gsr.ixSearched[1] + 1)*(gsr.iySearched[2] - gsr.iySearched[1] + 1)
end
