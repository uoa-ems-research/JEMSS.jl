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
	fnames = Set(fieldnames(Call))
	
	@assert(length(calls) == numCalls)
	@assert(length(backupCalls) == numCalls)
	
	# from fnames, remove fixed parameters
	fnamesFixed = Set([:index, :priority, :transport, :hospitalIndex, :location,
		:arrivalTime, :dispatchDelay, :onSceneDuration, :handoverDuration,
		:nearestNodeIndex, :nearestNodeDist])
	setdiff!(fnames, fnamesFixed)
	
	recentCallIndex = something(findlast(call -> call.arrivalTime <= sim.time, calls), 0)
	@assert(all(i -> calls[i].status == callNullStatus, recentCallIndex+1:numCalls))
	
	# reset calls that arrived before (or at) sim.time
	for fname in fnames
		if isprimitivetype(fieldtype(Call, fname))
			for i = 1:recentCallIndex
				setfield!(calls[i], fname, getfield(backupCalls[i], fname))
			end
		elseif fieldtype(Call, fname) == Location
			for i = 1:recentCallIndex
				copy!(getfield(calls[i], fname), getfield(backupCalls[i], fname)) # faster than deepcopy for Location
			end
		else
			for i = 1:recentCallIndex
				setfield!(calls[i], fname, deepcopy(getfield(backupCalls[i], fname)))
			end
		end
	end
end

function setCallStatus!(call::Call, status::CallStatus, time::Float)
	@assert(call.statusSetTime <= time)
	
	# stats for previous status
	statusDuration = time - call.statusSetTime
	prevStatus = call.status # shorthand
	if prevStatus == callQueued
		call.queuedDuration += statusDuration
	elseif prevStatus == callWaitingForAmb
		call.waitingForAmbDuration += statusDuration
		call.ambGoingToCallDuration = statusDuration # overwrite previous value if have attempted dispatch to call multiple times
	elseif prevStatus == callGoingToHospital
		@assert(call.transportDuration == 0.0) # should not have been calculated yet
		call.transportDuration = statusDuration
	end
	
	# stats for new status
	if status == callQueued
		call.wasQueued = true
	elseif status == callWaitingForAmb
		if prevStatus == callScreening
			@assert(isapprox(time, call.arrivalTime + call.dispatchDelay)) # if dispatching after call screening, only delay should be dispatch delay
		end
		call.dispatchTime = time
	elseif status == callOnSceneTreatment
		@assert(call.ambArrivalTime == nullTime) # value should not have been set yet
		call.ambArrivalTime = time
		call.responseDuration = time - call.arrivalTime
		call.bumpedDuration = call.waitingForAmbDuration - call.ambGoingToCallDuration
	elseif status == callAtHospital
		call.hospitalArrivalTime = time
	elseif status == callProcessed
		@assert(call.serviceDuration == 0.0) # value should not have been set yet
		call.serviceDuration = time - call.arrivalTime
	end
	
	call.status = status
	call.statusSetTime = time
end
