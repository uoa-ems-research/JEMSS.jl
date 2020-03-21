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

# Maximal Covering Location Problem (MCLP)

@testset "mclp" begin
	## tiny
	pointDemands = [1.0]
	coverMatrix = fill(true, 1, 1)
	results = Dict()
	@test solveMclp(0, pointDemands, coverMatrix, results = results) == [0]
	@test solveMclp(1, pointDemands, coverMatrix, results = results) == [1]
	
	## small
	pointDemands = Float[2, 2, 1, 1]
	coverMatrix = convert(Array{Bool,2}, [ # coverMatrix[i,j] = true if facility location i can cover point j
		1 1 0 0
		0 1 1 0
		1 0 0 1
	])
	
	@test solveMclp(0, pointDemands, coverMatrix, results = results) == [0, 0, 0]
	@test results[:objVal] == 0
	
	@test solveMclp(1, pointDemands, coverMatrix, results = results) == [1, 0, 0]
	@test results[:objVal] == 4
	
	@test solveMclp(2, pointDemands, coverMatrix, results = results) == [0, 1, 1]
	@test results[:objVal] == 6
	
	@test solveMclp(3, pointDemands, coverMatrix, results = results) == [1, 1, 1]
	@test results[:objVal] == 6
end
