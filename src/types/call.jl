# get next call from list
function getCall!(queuedCallList::Vector{Call})
	nextCall = nothing
	if length(queuedCallList) > 0
		nextCall = pop!(queuedCallList)
	end
	return nextCall
end

# add call to list that is sorted by call priorities and then by arrival times within priorities
function addCallToQueue!(queuedCallList::Vector{Call}, call::Call)
	# find where to insert new call into list
	# maintain sorting by priority (lower value means a higher priority),
	# and then sorting by arrival time within priority
	# calls nearer to end of list are higher priority and arrived sooner
	i = length(queuedCallList) + 1
	while i > 1 && (call.priority > queuedCallList[i-1].priority ||
		(call.priority == queuedCallList[i-1].priority && call.arrivalTime > queuedCallList[i-1].arrivalTime))
		i -= 1
	end
	insert!(queuedCallList, i, call)
	
	if checkMode
		# check that calls are ordered first by priority, then by time
		for i = length(queuedCallList): -1 : 2
			assert(queuedCallList[i].priority < queuedCallList[i-1].priority || 
				(queuedCallList[i].priority == queuedCallList[i-1].priority &&
				queuedCallList[i].arrivalTime <= queuedCallList[i-1].arrivalTime))
		end
	end
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
	
	# # from fnames, remove constant parameters
	# fnamesConst = Set([:index, :priority, :transfer, :location,
		# :arrivalTime, :dispatchDelay, :onSceneDuration, :transferDuration,
		# :nearestNodeIndex, :nearestNodeDist])
	# setdiff!(fnames, fnamesConst)
	
	for fname in fnames
		if typeof(getfield(nullCall, fname)) <: Number
			for i = 1:numCalls
				setfield!(calls[i], fname, getfield(backupCalls[i], fname))
			end
		else
			for i = 1:numCalls
				setfield!(calls[i], fname, deepcopy(getfield(backupCalls[i], fname)))
			end
		end
	end
end
