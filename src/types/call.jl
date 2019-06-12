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

# get next call from list
function getNextCall!(queuedCallList::Vector{Call})
	return length(queuedCallList) > 0 ? pop!(queuedCallList) : nothing
end

# reset simulation calls from sim.backup
# only reset calls with arrival time <= sim.time
# faster than sim.calls = deepcopy(sim.backup.calls)
function resetCalls!(sim::Simulation)
	@assert(!sim.backup.used)
	
	# shorthand:
	calls = sim.calls
	backupCalls = sim.backup.calls
	numCalls = sim.numCalls
	nullCall = Call()
	fnames = Set(fieldnames(Call))
	
	@assert(length(calls) == numCalls)
	@assert(length(backupCalls) == numCalls)
	
	# from fnames, remove fixed parameters
	fnamesFixed = Set([:index, :priority, :transport, :location,
		:arrivalTime, :dispatchDuration, :onSceneDuration, :handoverDuration,
		:nearestNodeIndex, :nearestNodeDist])
	setdiff!(fnames, fnamesFixed)
	
	recentCallIndex = something(findlast(call -> call.arrivalTime <= sim.time, calls), 0)
	@assert(all(i -> calls[i].status == callNullStatus, recentCallIndex+1:numCalls))
	
	# reset calls that arrived before (or at) sim.time
	for fname in fnames
		if typeof(getfield(nullCall, fname)) <: Number
			for i = 1:recentCallIndex
				setfield!(calls[i], fname, getfield(backupCalls[i], fname))
			end
		else
			for i = 1:recentCallIndex
				setfield!(calls[i], fname, deepcopy(getfield(backupCalls[i], fname)))
			end
		end
	end
end
