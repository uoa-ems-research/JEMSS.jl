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

# initialise given ambulance
# sets ambulance as sleeping, creates wake up event
function initAmbulance!(sim::Simulation, ambulance::Ambulance;
	wakeUpTime::Float = nullTime)
	wakeUpTime = (wakeUpTime == nullTime ? sim.startTime : wakeUpTime)
	
	@assert(ambulance.index != nullIndex)
	@assert(ambulance.stationIndex != nullIndex)
	@assert(sim.startTime <= wakeUpTime)
	
	ambulance.status = ambSleeping
	ambulance.prevStatusSetTime = sim.startTime
	# ambulance.stationIndex
	# ambulance.callIndex
	
	# initialise route to start at station
	ambulance.route = Route()
	station = sim.stations[ambulance.stationIndex]
	initRoute!(sim, ambulance.route; startLoc = station.location, startFNode = station.nearestNodeIndex, startFNodeDist = station.nearestNodeDist)
	
	# add wake up event
	addEvent!(sim.eventList; form = ambWakesUp, time = wakeUpTime, ambulance = ambulance)
end

# Set the ambulance status and time at which status started
# Mutates: ambulance
function setAmbStatus!(ambulance::Ambulance, status::AmbStatus, time::Float)
	# record duration of previous status
	@assert(ambulance.prevStatusSetTime <= time)
	ambulance.statusDurations[ambulance.status] += time - ambulance.prevStatusSetTime
	ambulance.prevStatusSetTime = time
	
	ambulance.status = status
end