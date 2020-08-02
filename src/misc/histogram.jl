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
	
	if (h1.closed == h2.closed) && (h1.isdensity == h2.isdensity) && (length(h1.edges) == 1 && length(h2.edges) == 1)
		e1, e2 = h1.edges[1], h2.edges[1] # shorthand
		if length(e1) < length(e2)
			h1, h2, e1, e2 = h2, h1, e2, e1 # swap h1 and h2, so that h1 is the 'longer' histogram
		end
		t1, t2 = typeof(e1), typeof(e2)
		# allow histograms to be added if they have the same bin size, with bins starting at the same point, but ending at different points
		if (t1 == t2) && (e1 != e2) && # note: if e1 == e2 here, then could just use merge(h1,h2)
			((t1 <: UnitRange && e1.start == e2.start) ||
			(t1 <: StepRange && e1.start == e2.start && e1.step == e2.step) ||
			(t1 <: StepRangeLen && e1.ref == e2.ref && e1.step == e2.step))
			
			h = deepcopy(h1)
			h.weights[1:length(h2.weights)] += h2.weights
			return h
		end
	end
	
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
