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

# using histogram from StatsBase

function Base.:+(h1::Histogram, h2::Histogram)::Histogram
	if h2 == nullHist return deepcopy(h1) end
	if h1 == nullHist return deepcopy(h2) end
	return merge(h1, h2) # requires that histograms have same field values, except for bin weights which just need to have equal length
end

function Base.:-(h::Histogram)::Histogram
	h2 = deepcopy(h)
	h2.weights = -h2.weights
	return h2
end

function Base.:-(h1::Histogram, h2::Histogram)::Histogram
	return h1 + (-h2)
end
