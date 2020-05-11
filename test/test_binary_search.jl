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

@testset "binary search" begin
	n = 10
	x = cumsum(rand(n).+0.1)
	@assert(issorted(x, lt=<=)) # values should be strictly increasing
	for i = 1:n
		@assert(binarySearch(x, x[i]*(1-eps())) == (i-1,i))
		@assert(binarySearch(x, x[i]) == (i,i))
		@assert(binarySearch(x, x[i]*(1+eps())) == (i,i+1))
	end
	@test true
end
