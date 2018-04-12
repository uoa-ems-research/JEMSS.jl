# get next call from list
function getCall!(queuedCallList::Vector{Call})
	nextCall = nothing
	if length(queuedCallList) > 0
		nextCall = pop!(queuedCallList)
	end
	return nextCall
end

# remove call from current calls, maintain sorting by call index
# return true if found, false otherwise
function removeCallFromCurrentCalls!(currentCallList::Vector{Call}, call::Call)
	# use linear search to find and then remove call
	# could use binary search since calls should be sorted by indices,
	# but there may not be enough calls in list to justify using it
	i = 1
	numCalls = length(currentCallList)
	while i <= numCalls
		if currentCallList[i].index == call.index
			# remove call from list
			deleteat!(currentCallList, i)
			return true
		end
		i += 1
	end
	
	return false
end

# reset simulation calls from sim.backup
# only reset calls with arrival time <= sim.time
# faster than sim.calls = deepcopy(sim.backup.calls)
function resetCalls!(sim::Simulation)
	assert(!sim.backup.used)
	
	# shorthand:
	calls = sim.calls
	backupCalls = sim.backup.calls
	numCalls = length(calls)
	nullCall = Call()
	fnames = Set(fieldnames(nullCall))
	
	assert(length(calls) == length(backupCalls))
	
	# from fnames, remove fixed parameters
	fnamesFixed = Set([:index, :priority, :transfer, :location,
		:arrivalTime, :dispatchDelay, :onSceneDuration, :transferDuration,
		:nearestNodeIndex, :nearestNodeDist])
	setdiff!(fnames, fnamesFixed)
	
	recentCallIndex = findlast(call -> call.arrivalTime <= sim.time, calls)
	assert(all(i -> calls[i].status == callNullStatus, recentCallIndex+1:numCalls))
	
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
