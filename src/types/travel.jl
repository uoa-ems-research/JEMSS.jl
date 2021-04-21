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

# get travel mode for given time and priority
# currentTime is sim.time, startTime is time that travel starts (>= sim.time)
function getTravelMode!(travel::Travel, priority::Priority, currentTime::Float; startTime::Float = currentTime)
	@assert(currentTime != nullTime)
	@assert(priority != nullPriority)
	@assert(currentTime <= startTime)
	
	travel.recentSetsStartTimesIndex = getTravelSetsStartTimesIndex(travel, currentTime)
	
	i = getTravelSetsStartTimesIndex(travel, startTime)
	travelSetIndex = travel.setsTimeOrder[i]
	travelModeIndex = travel.modeLookup[travelSetIndex, Int(priority)]
	
	return travel.modes[travelModeIndex]
end

# # get travel mode index for given time and priority
# function getTravelModeIndex!(travel::Travel, priority::Priority, startTime::Float)
	# travelMode = getTravelMode!(travel, priority, startTime)
	# return travelMode.index
# end

# get travel setsStartTimes index for given travel start time
function getTravelSetsStartTimesIndex(travel::Travel, startTime::Float)
	# shorthand:
	setsStartTimes = travel.setsStartTimes
	i = travel.recentSetsStartTimesIndex
	n = length(setsStartTimes)
	
	@assert(setsStartTimes[i] <= startTime) # otherwise, have gone back in time?
	
	# use linear search to find position in setsStartTimes
	while i < n && setsStartTimes[i+1] <= startTime
		i += 1
	end
	
	return i
end

# find nearest hospital to call location, given the travel priority
function nearestHospitalToCall!(sim::Simulation, call::Call, priority::Priority)
	travelMode = getTravelMode!(sim.travel, priority, sim.time)
	hospitalIndex = travelMode.fNetTravel.fNodeNearestHospitalIndex[call.nearestNodeIndex]
	return hospitalIndex
end
