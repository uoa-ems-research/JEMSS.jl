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
	fnames = Set(fieldnames(nullCall))
	
	@assert(length(calls) == numCalls)
	@assert(length(backupCalls) == numCalls)
	
	# from fnames, remove fixed parameters
	fnamesFixed = Set([:index, :priority, :transfer, :location,
		:arrivalTime, :dispatchDelay, :onSceneDuration, :transferDuration,
		:nearestNodeIndex, :nearestNodeDist])
	setdiff!(fnames, fnamesFixed)
	
	recentCallIndex = findlast(call -> call.arrivalTime <= sim.time, calls)
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
