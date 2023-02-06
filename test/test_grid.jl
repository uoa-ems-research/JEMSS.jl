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

using Random

@testset "grid" begin
    # create map
    mp = Map(0, 1, 0, 1, 10, 10) # Map(xMin, xMax, yMin, yMax, xScale, yScale)

    # create nodes
    numNodes = 100
    offRoadAccessProb = 0.8
    rng = Random.MersenneTwister(0)
    nodes = [Node() for i = 1:numNodes]
    for i = 1:numNodes
        nodes[i].index = i
        nodes[i].location = randLocation(mp; rng=rng)
        nodes[i].offRoadAccess = rand() < offRoadAccessProb
    end
    numOffRoadAccessNodes = sum(node -> node.offRoadAccess, nodes)
    @assert(numOffRoadAccessNodes < numNodes) # otherwise not all cases will be tested

    # create grids, add nodes to grid
    nx = ny = 10
    grid1 = Grid(mp, nx, ny) # grid with only nodes that have offRoadAccess == true
    JEMSS.gridPlaceNodes!(mp, grid1, nodes)
    grid2 = Grid(mp, nx, ny) # grid with all nodes
    JEMSS.gridPlaceNodes!(mp, grid2, nodes; offRoadAccessRequired=false)

    # test number of nodes added to each grid
    @test sum(gr -> length(gr.nodeIndices), grid1.rects) == numOffRoadAccessNodes
    @test sum(gr -> length(gr.nodeIndices), grid2.rects) == numNodes

    # each node should be closest to itself (for grid2, not grid1)
    @test all(node -> findNearestNode(mp, grid2, nodes, node.location)[1] == node.index, nodes)

    # check that grid rectangles have some different numbers of nodes, so that different cases will be tested
    @assert(all(grid -> any(rect -> length(rect.nodeIndices) == 0, grid.rects), (grid1, grid2)))
    @assert(all(grid -> any(rect -> length(rect.nodeIndices) == 1, grid.rects), (grid1, grid2)))
    @assert(all(grid -> any(rect -> length(rect.nodeIndices) >= 2, grid.rects), (grid1, grid2)))

    # compare findNearestNode function with brute-force method (no grid)
    grid1pass = grid2pass = true
    for x = range(0; length=10 * nx + 1, stop=1), y = range(0; length=10 * ny + 1, stop=1)
        location = Location(x, y)
        grid1pass &= findNearestNode(mp, grid1, nodes, location) == findNearestNode(mp, nodes, location; offRoadAccessRequired=true)
        grid2pass &= findNearestNode(mp, grid2, nodes, location) == findNearestNode(mp, nodes, location; offRoadAccessRequired=false)
    end
    @test grid1pass
    @test grid2pass

    # # plot of nearest node from sample locations
    # # this is not a test, but I wanted visual confirmation that the findNearestNode function was working
    # begin
    #     # sample locations
    #     xs = range(0; length=100 * nx + 1, stop=1)
    #     ys = range(0; length=100 * ny + 1, stop=1)
    #     z = zeros(Int, length(ys), length(xs))
    #     for (ix, x) in enumerate(xs), (iy, y) in enumerate(ys)
    #         location = Location(x, y)
    #         z[iy, ix] = findNearestNode(mp, grid2, nodes, location)[1]
    #     end

    #     using Plots
    #     Plots.plotly()
    #     plt = heatmap(xs, ys, z, aspect_ratio=1)
    #     x = [node.location.x for node in nodes]
    #     y = [node.location.y for node in nodes]
    #     Plots.scatter!(plt, x, y)
    # end
end
