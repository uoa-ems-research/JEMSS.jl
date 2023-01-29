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

# Run the function on the given items, using multi-threading.
# Kwarg `numThreads` can be set to reduce cpu usage, though this is only approximate.
# Requires environment variable JULIA_NUM_THREADS to be set before starting julia.
function runParallel!(f::Function, items...; numThreads::Int=maxNumThreads)
    @assert(numThreads >= 1)
    @threads for i = 1:numThreads
        for j = i:numThreads:length(items)
            f(items[j])
        end
    end
end
