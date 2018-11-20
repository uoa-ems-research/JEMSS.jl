# add call to list that is sorted by call priorities and then by arrival times within priorities
function addCallToQueueSortPriorityThenTime!(queuedCallList::Vector{Call}, call::Call)
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
			@assert(queuedCallList[i].priority < queuedCallList[i-1].priority ||
				(queuedCallList[i].priority == queuedCallList[i-1].priority &&
				queuedCallList[i].arrivalTime <= queuedCallList[i-1].arrivalTime))
		end
	end
end

