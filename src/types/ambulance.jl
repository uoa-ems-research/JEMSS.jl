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
	ambulance.statusSetTime = sim.startTime
	# ambulance.stationIndex
	# ambulance.callIndex
	
	# initialise route to start at station
	ambulance.route = Route()
	station = sim.stations[ambulance.stationIndex]
	initRoute!(sim, ambulance.route; startLoc = station.location, startFNode = station.nearestNodeIndex, startFNodeDist = station.nearestNodeDist)
	
	# statistics
	ambStatuses = setdiff(instances(AmbStatus), (ambNullStatus,))
	n = maximum(s -> Int(s), ambStatuses)
	@assert(all(s -> in(Int(s), 1:n), ambStatuses)) # statuses should be numbered 1:n
	ambulance.statusDurations = Dict([s => 0.0 for s in ambStatuses])
	ambulance.statusTransitionCounts = zeros(Int, n, n)
	
	# add wake up event
	addEvent!(sim.eventList; form = ambWakesUp, time = wakeUpTime, ambulance = ambulance)
end

# Set the ambulance status and time at which status started
# Mutates: ambulance
function setAmbStatus!(sim::Simulation, ambulance::Ambulance, status::AmbStatus, time::Float)
	@assert(ambulance.statusSetTime <= time)
	
	# stats
	statusDuration = time - ambulance.statusSetTime
	ambulance.statusDurations[ambulance.status] += statusDuration
	ambulance.statusTransitionCounts[Int(ambulance.status), Int(status)] += 1
	
	prevStatus = ambulance.status
	if isBusy(prevStatus)
		ambulance.totalBusyDuration += statusDuration
	end
	if isTravelling(prevStatus)
		@assert(time <= ambulance.route.endTime)
		ambulance.totalTravelDuration += time - ambulance.route.startTime
		ambulance.totalTravelDistance += calcRouteDistance!(sim, ambulance.route, time)
	end
	
	# stats
	if status == ambGoingToCall
		# dispatch stats
		ambulance.numDispatches += 1
		if prevStatus == ambIdleAtStation
			ambulance.numDispatchesFromStation += 1
		elseif isGoingToStation(prevStatus)
			ambulance.numDispatchesOnRoad += 1
		elseif prevStatus == ambGoingToCall
			ambulance.numDispatchesOnRoad += 1
			ambulance.numRedispatches += 1
		elseif prevStatus == ambFreeAfterCall
			ambulance.numDispatchesOnFree += 1
		else error()
		end
	elseif status == ambAtCall
		ambulance.numCallsTreated += 1
	elseif status == ambAtHospital
		ambulance.numCallsTransported += 1
	elseif status == ambMovingUpToStation
		# move up stats
		ambulance.numMoveUps += 1
		if prevStatus == ambIdleAtStation
			ambulance.numMoveUpsFromStation += 1
		elseif isGoingToStation(prevStatus)
			ambulance.numMoveUpsOnRoad += 1
		elseif prevStatus == ambFreeAfterCall
			ambulance.numMoveUpsOnFree += 1
		else error()
		end
	end
	
	ambulance.statusSetTime = time
	ambulance.prevStatus = ambulance.status
	ambulance.status = status
end

isBusy(s::AmbStatus)::Bool = in(s, (ambGoingToCall, ambAtCall, ambGoingToHospital, ambAtHospital))
isFree(s::AmbStatus)::Bool = in(s, (ambIdleAtStation, ambFreeAfterCall, ambReturningToStation, ambMovingUpToStation))
isWorking(s::AmbStatus)::Bool = !in(s, (ambNullStatus, ambSleeping)) # should return same value as isBusy(s) || isFree(s)
isGoingToStation(s::AmbStatus)::Bool = in(s, (ambReturningToStation, ambMovingUpToStation))
isTravelling(s::AmbStatus)::Bool = in(s, (ambGoingToCall, ambGoingToHospital, ambReturningToStation, ambMovingUpToStation))
