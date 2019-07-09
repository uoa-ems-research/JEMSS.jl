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

@info("Todo: check correctness of statistics.")

# Capture statistics at current time.
# Mutates: sim.stats
function captureSimStats!(sim::Simulation, currentTime::Float)
	if sim.complete currentTime = min(currentTime, sim.endTime) end
	
	if !sim.complete @assert(currentTime == sim.stats.nextCaptureTime) end
	@assert(sim.time <= currentTime)
	@assert(isempty(sim.eventList) || currentTime <= sim.eventList[end].time) # currentTime should be before next event
	
	capture = SimPeriodStats()
	capture.startTime = sim.startTime
	capture.endTime = currentTime
	capture.duration = capture.endTime - capture.startTime
	capture.ambulances = [AmbulanceStats(sim, a) for a in sim.ambulances]
	capture.hospitals = [HospitalStats(h) for h in sim.hospitals]
	# capture.stations = [StationStats(s) for s in sim.stations]
	capture.ambulance = sum(capture.ambulances)
	capture.hospital = sum(capture.hospitals)
	# capture.station = sum(capture.stations)
	
	push!(sim.stats.captures, capture)
	
	# determine next capture time
	i = sim.stats.nextCaptureTimeIndex += 1
	sim.stats.nextCaptureTime = i <= length(sim.stats.captureTimes) ? sim.stats.captureTimes[i] : Inf
end

# For statistics that should be calculated at sim end.
# Mutates: sim.stats
function populateSimStats!(sim::Simulation)
	@assert(sim.complete)
	
	# shorthand
	stats = sim.stats
	captures = stats.captures
	
	@assert(all(s -> s.startTime == sim.startTime, captures))
	@assert(issorted(captures, by = s -> s.endTime))
	@assert(captures[end].endTime == sim.endTime)
	
	# populate call statistics
	callsByPriority = Dict([p => filter(c -> c.priority == p, sim.calls) for p in priorities])
	for i = 1:length(captures)
		period = i == 1 ? captures[i] : captures[i] - captures[i-1]
		period.call = CallStats(sim, sim.calls, period.startTime, period.endTime)
		for priority in priorities
			period.callPriorities[priority] = CallStats(sim, callsByPriority[priority], period.startTime, period.endTime)
		end
		push!(stats.periods, period)
		
		captures[i].call = i == 1 ? period.call : period.call + captures[i-1].call
		captures[i].callPriorities = i == 1 ? period.callPriorities : period.callPriorities + captures[i-1].callPriorities
	end
end

# Return ambulance stats for given ambulance.
# Accounts for durations and distances of partially completed travelling and processes.
function AmbulanceStats(sim::Simulation, ambulance::Ambulance)::AmbulanceStats
	stats = AmbulanceStats()
	
	# copy some fields
	stats.ambIndex = ambulance.index
	stats.numCallsTreated = ambulance.numCallsTreated
	stats.numCallsTransported = ambulance.numCallsTransported
	stats.numDispatchesFromStation = ambulance.numDispatchesFromStation
	stats.numDispatchesOnRoad = ambulance.numDispatchesOnRoad
	stats.numDispatchesOnFree = ambulance.numDispatchesOnFree
	stats.numRedispatches = ambulance.numRedispatches
	stats.numMoveUpsFromStation = ambulance.numMoveUpsFromStation
	stats.numMoveUpsOnRoad = ambulance.numMoveUpsOnRoad
	stats.numMoveUpsOnFree = ambulance.numMoveUpsOnFree
	stats.statusTransitionCounts = deepcopy(ambulance.statusTransitionCounts)
	
	stats.statusDurations = deepcopy(ambulance.statusDurations)
	stats.statusDurations[ambulance.status] += sim.time - ambulance.statusSetTime # to make sure that sum of statusDurations equals the sim duration
	
	# calculate
	stats.numDispatches = stats.numDispatchesFromStation + stats.numDispatchesOnRoad + stats.numDispatchesOnFree + stats.numRedispatches
	stats.numMoveUps = stats.numMoveUpsFromStation + stats.numMoveUpsOnRoad + stats.numMoveUpsOnFree
	
	# calculate travel distance, accounting for any partially completed route
	stats.totalTravelDistance = ambulance.totalTravelDistance # ambulance.totalTravelDistance is for routes completed before sim.time
	if sim.time < ambulance.route.endTime # have not finished route
		stats.totalTravelDistance += calcRouteDistance!(sim, ambulance.route, sim.time)
	end
	
	busyStatuses = (ambGoingToCall, ambAtCall, ambGoingToHospital, ambAtHospital)
	travelStatuses = (ambGoingToCall, ambGoingToHospital, ambGoingToStation, ambMovingUpToStation)
	stats.totalBusyDuration = sum(status -> stats.statusDurations[status], busyStatuses) # >= ambulance.totalBusyDuration
	stats.totalTravelDuration = sum(status -> stats.statusDurations[status], travelStatuses) # >= ambulance.totalTravelDuration
	
	if checkMode
		@assert(isapprox(sum(values(stats.statusDurations)), sim.time - sim.startTime))
		@assert(stats.totalBusyDuration >= ambulance.totalBusyDuration || isapprox(stats.totalBusyDuration, ambulance.totalBusyDuration))
		@assert(stats.totalTravelDuration >= ambulance.totalTravelDuration || isapprox(stats.totalTravelDuration, ambulance.totalTravelDuration))
	end
	
	return stats
end

function CallStats(sim::Simulation, calls::Vector{Call})::CallStats
	if isempty(calls) return CallStats() end
	@assert(calls[end].status == callProcessed) # should be true for preceeding calls also
	stats = CallStats()
	stats.numCalls = length(calls)
	if stats.numCalls == 1 stats.callIndex = calls[1].index end
	stats.numQueued = count(c -> c.arrivalTime + c.dispatchDelay < c.dispatchTime, calls) # assuming no mobilisation delay
	stats.numBumped = count(c -> c.numBumps > 0, calls)
	stats.numBumps = sum(c -> c.numBumps, calls)
	stats.numTransports = count(c -> c.transport, calls)
	stats.numResponsesInTime = count(c -> c.responseDuration <= sim.targetResponseDurations[Int(c.priority)], calls)
	stats.totalResponseDuration = sum(c -> c.responseDuration, calls)
	stats.totalQueuedDuration = sum(c -> c.queuedDuration, calls)
	stats.totalOnSceneDuration = sum(c -> c.onSceneDuration, calls)
	stats.totalTransportDuration = sum(c -> c.transport ? c.hospitalArrivalTime - (c.ambArrivalTime + c.onSceneDuration) : 0, calls)
	stats.totalAtHospitalDuration = sum(c -> c.transport ? c.handoverDuration : 0, calls)
	# stats.totalDispatchDelay = sum(c -> c.dispatchTime - c.arrivalTime, calls) # can be different from dispatchDelay due to queueing and ambulance redispatching
	return stats
end

function CallStats(sim::Simulation, calls::Vector{Call}, startTime::Float, endTime::Float)::CallStats
	# filter calls by arrivalTime within [startTime, endTime)
	i = findfirst(c -> c.arrivalTime >= startTime, calls)
	if i == nothing return CallStats() end # no calls in time period
	j = something(findnext(c -> c.arrivalTime >= endTime, calls, i), length(calls) + 1) - 1
	return CallStats(sim, calls[i:j])
end

function HospitalStats(hospital::Hospital)::HospitalStats
	stats = HospitalStats()
	stats.hospitalIndex = hospital.index
	stats.numCalls = hospital.numCalls
	return stats
end

function StationStats(station::Station)::StationStats
	stats = StationStats()
	stats.stationIndex = station.index
	# todo
	return stats
end

# for ambulance statusDurations:
Base.:+(a::Dict{AmbStatus,Float}, b::Dict{AmbStatus,Float}) = merge(+, a, b)
Base.:-(a::Dict{AmbStatus,Float}, b::Dict{AmbStatus,Float}) = merge(-, a, b)

function Base.:+(a::AmbulanceStats, b::AmbulanceStats)::AmbulanceStats
	stats = AmbulanceStats()
	for fname in fieldnames(AmbulanceStats)
		setfield!(stats, fname, getfield(a, fname) + getfield(b, fname))
	end
	stats.ambIndex = a.ambIndex == b.ambIndex ? a.ambIndex : nullIndex
	return stats
end

function Base.:-(a::AmbulanceStats, b::AmbulanceStats)::AmbulanceStats
	stats = AmbulanceStats()
	for fname in fieldnames(AmbulanceStats)
		setfield!(stats, fname, getfield(a, fname) - getfield(b, fname))
	end
	stats.ambIndex = a.ambIndex == b.ambIndex ? a.ambIndex : nullIndex
	return stats
end

function Base.:+(a::CallStats, b::CallStats)::CallStats
	stats = CallStats()
	for fname in fieldnames(CallStats)
		setfield!(stats, fname, getfield(a, fname) + getfield(b, fname))
	end
	stats.callIndex = a.callIndex == b.callIndex ? a.callIndex : nullIndex
	return stats
end

function Base.:-(a::CallStats, b::CallStats)::CallStats
	stats = CallStats()
	for fname in fieldnames(CallStats)
		setfield!(stats, fname, getfield(a, fname) - getfield(b, fname))
	end
	stats.callIndex = a.callIndex == b.callIndex ? a.callIndex : nullIndex
	return stats
end

Base.:+(a::Dict{Priority,CallStats}, b::Dict{Priority,CallStats}) = merge(+, a, b)
Base.:-(a::Dict{Priority,CallStats}, b::Dict{Priority,CallStats}) = merge(-, a, b)

function Base.:+(a::HospitalStats, b::HospitalStats)::HospitalStats
	stats = HospitalStats()
	for fname in fieldnames(HospitalStats)
		setfield!(stats, fname, getfield(a, fname) + getfield(b, fname))
	end
	stats.hospitalIndex = a.hospitalIndex == b.hospitalIndex ? a.hospitalIndex : nullIndex
	return stats
end

function Base.:-(a::HospitalStats, b::HospitalStats)::HospitalStats
	stats = HospitalStats()
	for fname in fieldnames(HospitalStats)
		setfield!(stats, fname, getfield(a, fname) - getfield(b, fname))
	end
	stats.hospitalIndex = a.hospitalIndex == b.hospitalIndex ? a.hospitalIndex : nullIndex
	return stats
end

function Base.:+(a::StationStats, b::StationStats)::StationStats
	stats = StationStats()
	for fname in fieldnames(StationStats)
		setfield!(stats, fname, getfield(a, fname) + getfield(b, fname))
	end
	stats.stationIndex = a.stationIndex == b.stationIndex ? a.stationIndex : nullIndex
	return stats
end

function Base.:-(a::StationStats, b::StationStats)::StationStats
	stats = StationStats()
	for fname in fieldnames(StationStats)
		setfield!(stats, fname, getfield(a, fname) - getfield(b, fname))
	end
	stats.stationIndex = a.stationIndex == b.stationIndex ? a.stationIndex : nullIndex
	return stats
end

function Base.:+(a::SimPeriodStats, b::SimPeriodStats)::SimPeriodStats
	@assert(a.startTime == b.endTime || a.endTime == b.startTime)
	stats = SimPeriodStats()
	for fname in fieldnames(SimPeriodStats)
		setfield!(stats, fname, getfield(a, fname) + getfield(b, fname))
	end
	stats.startTime = min(a.startTime, b.startTime)
	stats.endTime = max(a.endTime, b.endTime)
	return stats
end

function Base.:-(a::SimPeriodStats, b::SimPeriodStats)::SimPeriodStats
	@assert(a.startTime == b.startTime)
	stats = SimPeriodStats()
	for fname in fieldnames(SimPeriodStats)
		setfield!(stats, fname, getfield(a, fname) - getfield(b, fname))
	end
	stats.startTime = b.endTime
	stats.endTime = a.endTime
	return stats
end
