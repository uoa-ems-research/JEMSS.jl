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

@testset "histogram addition" begin
	# null histogram
	h1 = fit(Histogram, [0], 0:1)
	h2 = NullHist()
	h3 = h1+h2
	@test h3 == h1
	@test h3 == h2+h1 # commutative
	
	# histograms with UnitRange edges
	h1 = fit(Histogram, [0], 0:2)
	h2 = fit(Histogram, [0,2], 0:3)
	h3 = h1+h2
	@test h3 == h2+h1 # commutative
	@test h3.edges == h2.edges
	@test h3.weights == [2,0,1] # [1,0] + [1,0,1]
	
	# histograms with StepRange edges
	h1 = fit(Histogram, [0], 0:2:4)
	h2 = fit(Histogram, [0,5], 0:2:6)
	h3 = h1+h2
	@test h3 == h2+h1 # commutative
	@test h3.edges == h2.edges
	@test h3.weights == [2,0,1] # [1,0] + [1,0,1]
	
	# histograms with StepRangeLen edges
	h1 = fit(Histogram, [0], 0.0:2.0)
	h2 = fit(Histogram, [0,2], 0.0:3.0)
	h3 = h1+h2
	@test h3 == h2+h1 # commutative
	@test h3.edges == h2.edges
	@test h3.weights == [2,0,1] # [1,0] + [1,0,1]
end
