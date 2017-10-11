# get travel mode for given time and priority
function getTravelMode!(travel::Travel, priority::Priority, startTime::Float)
	assert(startTime != nullTime)
	assert(priority != nullPriority)
	
	updateTravelToTime!(travel, startTime) # update travel.recentSetsStartTimesIndex for given start time
	travelSetIndex = travel.setsTimeOrder[travel.recentSetsStartTimesIndex]
	travelModeIndex = travel.modeLookup[travelSetIndex, Int(priority)]
	
	return travel.modes[travelModeIndex]
end

# # get travel mode index for given time and priority
# function getTravelModeIndex!(travel::Travel, priority::Priority, startTime::Float)
	# travelMode = getTravelMode!(travel, priority, startTime)
	# return travelMode.index
# end

# update travel.recentSetsStartTimesIndex up to given travel start time
function updateTravelToTime!(travel::Travel, startTime::Float)
	# shorthand:
	setsStartTimes = travel.setsStartTimes
	i = travel.recentSetsStartTimesIndex
	n = length(setsStartTimes)
	
	assert(setsStartTimes[i] <= startTime) # otherwise, have gone back in time?
	
	if i == n || startTime < setsStartTimes[i+1]
		return # do nothing
	end
	
	# use linear search to update position in setsStartTimes
	while i < n && setsStartTimes[i+1] <= startTime
		i += 1
	end
	
	travel.recentSetsStartTimesIndex = i
end
