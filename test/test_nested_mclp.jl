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

# Maximal Coverage Location Problem applied to creating a nested compliance table

@testset "nested mclp" begin
    ## tiny
    pointDemands = [1.0]
    coverMatrix = fill(true, 1, 1)
    weights = [1.0]
    results = Dict()
    @test solveNestedMclp(pointDemands, coverMatrix, weights; results=results) == [1]

    ## small
    pointDemands = Float[5, 3, 2, 1]
    coverMatrix = convert(
        Array{Bool,2},
        [ # coverMatrix[i,j] = true if facility location i can cover point j
            1 1 0 0
            0 1 1 0
            1 0 0 1
        ]
    )

    weights = [1, 1, 1] / 3
    @test solveNestedMclp(pointDemands, coverMatrix, weights, results=results) == [1, 2, 3]
    results[:objVal] = sum(weights .* [5 + 3, 5 + 3 + 2, 5 + 3 + 2 + 1])

    weights = [1, 4, 1] / 6
    @test solveNestedMclp(pointDemands, coverMatrix, weights, results=results) == [3, 2, 1]
    results[:objVal] = sum(weights .* [5 + 1, 5 + 3 + 2 + 1, 5 + 3 + 2 + 1])
end
