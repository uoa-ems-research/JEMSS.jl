# Demand (call) model

# get demand mode for given time and priority
function getDemandMode!(demand::Demand, priority::Priority, startTime::Float)
	@assert(startTime != nullTime)
	@assert(priority != nullPriority)
	
	updateDemandToTime!(demand, startTime) # update demand.recentSetsStartTimesIndex for given start time
	demandSetIndex = demand.setsTimeOrder[demand.recentSetsStartTimesIndex]
	demandModeIndex = demand.modeLookup[demandSetIndex, Int(priority)]
	
	return demand.modes[demandModeIndex]
end

# update demand.recentSetsStartTimesIndex up to given demand start time
function updateDemandToTime!(demand::Demand, startTime::Float)
	# shorthand:
	setsStartTimes = demand.setsStartTimes
	i = demand.recentSetsStartTimesIndex
	n = length(setsStartTimes)
	
	@assert(setsStartTimes[i] <= startTime) # otherwise, have gone back in time?
	
	if i == n || startTime < setsStartTimes[i+1]
		return # do nothing
	end
	
	# use linear search to update position in setsStartTimes
	while i < n && setsStartTimes[i+1] <= startTime
		i += 1
	end
	
	demand.recentSetsStartTimesIndex = i
end
