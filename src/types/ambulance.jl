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
	
	# initialise route to start at station
	ambulance.route = Route()
	station = sim.stations[ambulance.stationIndex]
	initRoute!(sim, ambulance.route; startLoc = station.location, startFNode = station.nearestNodeIndex, startFNodeDist = station.nearestNodeDist)
	
	# statistics
	ambStatuses = setdiff(instances(AmbStatus), (ambNullStatus,))
	n = maximum(s -> Int(s), ambStatuses)
	@assert(all(s -> in(Int(s), 1:n), ambStatuses)) # statuses should be numbered 1:n
	ambulance.statusTransitionCounts = zeros(Int, n, n)
	statusesAndSets = (ambStatuses..., instances(AmbStatusSet)...)
	ambulance.statusDurations = Dict([s => 0.0 for s in statusesAndSets])
	ambulance.statusDistances = Dict([s => 0.0 for s in (ambStatusSets[ambTravelling]..., instances(AmbStatusSet)...)])
	
	# add wake up event
	addEvent!(sim.eventList; form = ambWakesUp, time = wakeUpTime, ambulance = ambulance, station = station)
end

function addStatusDuration!(statusDurations::Dict{Union{AmbStatus,AmbStatusSet},Float}, status::AmbStatus, duration::Float)
	statusDurations[status] += duration
	for s in ambStatusToSets[status] statusDurations[s] += duration end # for status sets
end

function addStatusDistance!(statusDistances::Dict{Union{AmbStatus,AmbStatusSet},Float}, status::AmbStatus, dist::Float)
	@assert(isTravelling(status))
	statusDistances[status] += dist
	for s in ambStatusToSets[status] statusDistances[s] += dist end # for status sets
end

# Set the ambulance status and time at which status started
# Mutates: ambulance
function setAmbStatus!(sim::Simulation, ambulance::Ambulance, status::AmbStatus, time::Float)
	@assert(ambulance.statusSetTime <= time)
	
	# stats - previous status duration
	prevStatus = ambulance.status
	statusDuration = time - ambulance.statusSetTime
	addStatusDuration!(ambulance.statusDurations, prevStatus, statusDuration)
	if isTravelling(prevStatus)
		@assert(ambulance.statusSetTime == ambulance.route.startTime)
		@assert(time <= ambulance.route.endTime) # status should change when route ends, not after
		dist = calcRouteDistance!(sim, ambulance.route, time)
		addStatusDistance!(ambulance.statusDistances, prevStatus, dist)
	end
	
	ambulance.statusSetTime = time
	ambulance.prevStatus = prevStatus
	ambulance.status = status
	
	if sim.complete
		@assert(time == sim.endTime)
		@assert(prevStatus == status)
		return
	end
	
	ambulance.statusTransitionCounts[Int(prevStatus), Int(status)] += 1
	
	# stats
	if status == ambMobilising
		# dispatch stats
		ambulance.numDispatches += 1
		if prevStatus == ambIdleAtStation
			ambulance.numDispatchesFromStation += 1
		elseif prevStatus == ambMobilising # can happen for redispatch
			ambulance.numDispatchesWhileMobilising += 1
			ambulance.numRedispatches += 1
		else error()
		end
	elseif status == ambGoingToCall
		# dispatch stats
		if prevStatus != ambMobilising # already counted dispatches where prevStatus == ambMobilising
			ambulance.numDispatches += 1
		end
		if prevStatus == ambIdleAtStation
			ambulance.numDispatchesFromStation += 1
		elseif prevStatus == ambMobilising
		# 	ambulance.numDispatchesWhileMobilising += 1 # removed; this is already counted
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
	
	if status != ambMovingUpToStation
		ambulance.moveUpFromStationIndex = nullIndex
	end
end

isBusy(s::AmbStatus)::Bool = in(s, (ambMobilising, ambGoingToCall, ambAtCall, ambGoingToHospital, ambAtHospital))
isFree(s::AmbStatus)::Bool = in(s, (ambIdleAtStation, ambFreeAfterCall, ambReturningToStation, ambMovingUpToStation))
isWorking(s::AmbStatus)::Bool = !in(s, (ambNullStatus, ambSleeping)) # ambulance is operational, whether busy or free; should return same value as isBusy(s) || isFree(s)
isGoingToStation(s::AmbStatus)::Bool = in(s, (ambReturningToStation, ambMovingUpToStation))
isTravelling(s::AmbStatus)::Bool = in(s, (ambGoingToCall, ambGoingToHospital)) || isGoingToStation(s)

ambStatuses = collect(instances(AmbStatus))
const ambStatusSets = Dict{AmbStatusSet,Set{AmbStatus}}(
	ambWorking => Set(filter(s -> isWorking(s), ambStatuses)),
	ambBusy => Set(filter(s -> isBusy(s), ambStatuses)),
	ambFree => Set(filter(s -> isFree(s), ambStatuses)),
	ambTravelling => Set(filter(s -> isTravelling(s), ambStatuses)),
	ambGoingToStation => Set(filter(s -> isGoingToStation(s), ambStatuses))
)

ambStatusToSets = Dict([status => Set{AmbStatusSet}() for status in instances(AmbStatus)])
for (k, v) in ambStatusSets, ambStatus in v
	push!(ambStatusToSets[ambStatus], k)
end
