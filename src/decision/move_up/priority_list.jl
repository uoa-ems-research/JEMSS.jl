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

# priority list, for newly freed ambulances only
# (a priority list which uses all free ambulances can be described as a compliance table)

# initialise data relevant to move up
function initPriorityList!(sim::Simulation, priorityListFilename::String)
	# read priority list
	priorityList = readPriorityListFile(priorityListFilename)
	initPriorityList!(sim, priorityList)
end

# initialise data relevant to move up
function initPriorityList!(sim::Simulation, priorityList::Vector{Int})
	# shorthand names:
	ambulances = sim.ambulances
	stations = sim.stations
	numAmbs = sim.numAmbs
	numStations = sim.numStations
	pld = sim.moveUpData.priorityListData
	
	pld.priorityList = priorityList
	@assert(length(pld.priorityList) == numAmbs)
	@assert(maximum(pld.priorityList) <= numStations)
	
	# check that priority list does not violate station capacity constraints
	# @assert(all( sum(pld.priorityList .== j) .<= stations[j].capacity for j = 1:numStations ))
	for j = 1:numStations
		@assert(sum(pld.priorityList .== j) <= stations[j].capacity)
	end
	
	pld.stationNumIdleAmbs = Vector{Int}(undef, numStations)
	
	# check that ambulance to station assignments follow priority list
	pld.stationNumIdleAmbs[:] .= 0
	for i = 1:numAmbs
		pld.stationNumIdleAmbs[ambulances[i].stationIndex] += 1
	end
	if !all(pld.stationNumIdleAmbs[i] == sum(pld.priorityList .== i) for i = 1:numStations)
		@warn("Number of ambulances at each station does not match priority list")
	end
end

function priorityListMoveUp(sim::Simulation, newlyIdleAmb::Ambulance)
	@assert(sim.moveUpData.useMoveUp)
	
	# shorthand:
	pld = sim.moveUpData.priorityListData
	priorityList = pld.priorityList
	stationNumIdleAmbs = pld.stationNumIdleAmbs
	ambulances = sim.ambulances
	stations = sim.stations
	numAmbs = sim.numAmbs
	
	# calculate the number of idle ambulances at (or travelling to) each station
	stationNumIdleAmbs[:] .= 0
	for i = 1:numAmbs
		# do not count newly idle ambulance, it has not been assigned a station
		if isAmbAvailableForMoveUp(ambulances[i]) && i != newlyIdleAmb.index
			stationNumIdleAmbs[ambulances[i].stationIndex] += 1
		end
	end
	numIdleAmbs = sum(stationNumIdleAmbs) + 1
	
	# return first item in priority list that is not already filled
	stationIndex = nullIndex
	for i = 1:numIdleAmbs
		j = priorityList[i] # station index
		stationNumIdleAmbs[j] -= 1
		if stationNumIdleAmbs[j] < 0
			stationIndex = j
			break
		end
	end
	@assert(stationIndex != nullIndex)
	
	return [newlyIdleAmb], [stations[stationIndex]]
end
