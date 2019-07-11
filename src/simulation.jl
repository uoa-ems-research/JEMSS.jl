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

# run simulation

"""
	function simulate!(sim::Simulation;
		time::Float = Inf, duration::Float = Inf, numEvents::Int = -1,
		doPrint::Bool = false, printingInterval::Float = 1.0)
Run the simulation to completion, or to any specified stopping point (according to `time`, `duration`, or `numEvents`), whichever comes first. Returns `true` if simulation is complete; `false` otherwise.

# Keyword arguments
- `time` is the maximum time to simulate to; will stop if time of next event is greater
- `duration` is the maximum duration to simulate for; equivalent to setting `time = sim.time + duration`
- `numEvents` is the maximum number of events to simulate
- `doPrint` controls whether any lines are printed while simulating
- `printingInterval` is the time interval (seconds) between printing while simulating, if `doPrint = true`
"""
function simulate!(sim::Simulation;
	time::Real = Inf, duration::Real = Inf, numEvents::Real = Inf,
	doPrint::Bool = false, printingInterval::Real = 1.0)
	
	@assert(time == Inf || duration == Inf, "can only set one of: time, duration")
	@assert(time >= sim.time)
	@assert(duration >= 0)
	@assert(numEvents >= 0)
	@assert(printingInterval >= 1e-2, "printing too often can slow the simulation")
	
	duration != Inf && (time = sim.time + duration) # set end time based on duration
	
	# for printing progress
	startTime = Base.time()
	doPrint || (printingInterval = Inf)
	nextPrintTime = startTime + printingInterval
	eventCount = 0
	printProgress() = doPrint && print(@sprintf("\rsim time: %-9.2f sim duration: %-9.2f events simulated: %-9d real duration: %.2f seconds", sim.time, sim.time - sim.startTime, eventCount, Base.time() - startTime))
	
	# simulate
	stats = sim.stats # shorthand
	while !sim.complete && sim.eventList[end].time <= time && eventCount < numEvents
		if stats.doCapture && stats.nextCaptureTime <= sim.eventList[end].time
			captureSimStats!(sim, stats.nextCaptureTime)
		end
		
		simulateNextEvent!(sim)
		eventCount += 1
		
		if doPrint && Base.time() >= nextPrintTime
			printProgress()
			nextPrintTime += printingInterval
		end
	end
	
	printProgress()
	doPrint && println()
	
	if stats.doCapture && sim.complete
		captureSimStats!(sim, sim.endTime)
		populateSimStats!(sim)
	end
	
	return sim.complete
end

# simulate up to last event with time <= given time
# similar to calling simulate!(sim, time = time), but does not record statistics
# kept this for animation, to avoid problems with animating as an asynchronous task
function simulateToTime!(sim::Simulation, time::Float)
	while !sim.complete && sim.eventList[end].time <= time # eventList must remain sorted by non-increasing time
		simulateNextEvent!(sim)
	end
end

# run simulation to end
# similar to calling simulate!(sim), but does not record statistics
function simulateToEnd!(sim::Simulation)
	while !sim.complete
		simulateNextEvent!(sim)
	end
end

# set sim.backup to copy of sim (before running)
# note that for reducing memory usage, sim.backup does not contain backups of all fields
function backup!(sim::Simulation)
	@assert(!sim.used)
	
	# remove certain fields from sim before copying sim
	(net, travel, grid, resim, demand, demandCoverage) = (sim.net, sim.travel, sim.grid, sim.resim, sim.demand, sim.demandCoverage)
	(sim.net, sim.travel, sim.grid, sim.resim, sim.demand, sim.demandCoverage) = (Network(), Travel(), Grid(), Resimulation(), Demand(), DemandCoverage())
	
	sim.backup = deepcopy(sim)
	
	(sim.net, sim.travel, sim.grid, sim.resim, sim.demand, sim.demandCoverage) = (net, travel, grid, resim, demand, demandCoverage)
end
backupSim!(sim::Simulation) = backup!(sim) # compat

# reset sim from sim.backup
function reset!(sim::Simulation)
	@assert(!sim.backup.used)
	
	if sim.used
		resetCalls!(sim) # reset calls from sim.backup, need to do this before resetting sim.time
		
		fnames = Set(fieldnames(Simulation))
		fnamesDontCopy = Set([:backup, :net, :travel, :grid, :resim, :calls, :demand, :demandCoverage, :animating]) # will not (yet) copy these fields from sim.backup to sim
		# note that sim.backup does not contain a backup of all fields
		setdiff!(fnames, fnamesDontCopy) # remove fnamesDontCopy from fnames
		for fname in fnames
			try
				setfield!(sim, fname, deepcopy(getfield(sim.backup, fname)))
			catch e
				error(e)
			end
		end
		
		# reset resimulation state
		sim.resim.prevEventIndex = 0
		
		# reset travel and demand state
		sim.travel.recentSetsStartTimesIndex = 1
		sim.demand.recentSetsStartTimesIndex = 1
	end
end
resetSim!(sim::Simulation) = reset!(sim) # compat

# simulate next event in list
function simulateNextEvent!(sim::Simulation)
	# get next event, update event index and sim time
	event = getNextEvent!(sim.eventList)
	if event.form == nullEvent
		error()
	end
	sim.used = true
	sim.eventIndex += 1
	event.index = sim.eventIndex
	sim.time = event.time
	
	if sim.resim.use
		resimCheckCurrentEvent!(sim, event)
	elseif sim.writeOutput
		writeEventToFile!(sim, event)
	end
	
	simulateEvent!(sim, event)
	
	if length(sim.eventList) == 0
		# simulation complete
		@assert(sim.endTime == nullTime)
		@assert(sim.complete == false)
		sim.endTime = sim.time
		sim.complete = true
	end
end

# simulate event
function simulateEvent!(sim::Simulation, event::Event)
	# format:
	# next event may change relevant ambulance / call fields at event.time
	# event may then trigger future events / cancel scheduled events
	
	@assert(sim.time == event.time)
	
	eventForm = event.form
	if eventForm == nullEvent
		error("null event")
		
################
	
	elseif eventForm == ambGoesToSleep
		ambulance = sim.ambulances[event.ambIndex]
		@assert(ambulance.status == ambIdleAtStation) # ambulance should be at station before sleeping
		@assert(event.callIndex == nullIndex)
		
		setAmbStatus!(ambulance, ambSleeping, sim.time)
		# ambulance.stationIndex
		# ambulance.callIndex
		
		addEvent!(sim.eventList; parentEvent = event, form = ambWakesUp, time = sim.time + sleepDuration, ambulance = ambulance)
		
################
	
	elseif eventForm == ambWakesUp
		ambulance = sim.ambulances[event.ambIndex]
		@assert(ambulance.status == ambSleeping) # ambulance should have been sleeping
		@assert(event.callIndex == nullIndex)
		
		setAmbStatus!(ambulance, ambIdleAtStation, sim.time)
		# ambulance.stationIndex
		# ambulance.callIndex
		
################
	
	elseif eventForm == callArrives
		@assert(event.ambIndex == nullIndex) # no ambulance should be assigned yet
		call = sim.calls[event.callIndex]
		@assert(call.status == callNullStatus)
		
		push!(sim.currentCalls, sim.calls[event.callIndex])
		call.status = callScreening
		
		addEvent!(sim.eventList; parentEvent = event, form = considerDispatch, time = sim.time + call.dispatchDelay, call = call)
		
		# add next call arrival to event queue
		if call.index < sim.numCalls
			addEvent!(sim.eventList, sim.calls[call.index + 1])
		end
		
################
	
	elseif eventForm == considerDispatch
		@assert(event.ambIndex == nullIndex) # no ambulance should be assigned yet
		call = sim.calls[event.callIndex]
		@assert(call.status == callScreening || call.status == callWaitingForAmb) # callWaitingForAmb if call bumped
		
		# find ambulance to respond
		if sim.resim.use
			ambIndex = resimFindAmbToDispatch(sim, call)
		else
			ambIndex = sim.findAmbToDispatch!(sim, call)
		end
		# if no ambulance found, add call to queue,
		# will respond when an ambulance becomes free
		if ambIndex == nullIndex
			sim.addCallToQueue!(sim.queuedCallList, call)
			call.status = callQueued
			call.wasQueued = true # stats
		else
			ambulance = sim.ambulances[ambIndex]
			
			# remove any previously scheduled events
			ambWasOnRoute = true
			if isGoingToStation(ambulance.status)
				deleteEvent!(sim.eventList, ambulance.event)
			elseif ambulance.status == ambGoingToCall
				# if ambulance was going to call, redirect ambulance
				deleteEvent!(sim.eventList, ambulance.event)
				
				# bumped call has no ambulance assigned yet
				sim.calls[ambulance.callIndex].ambIndex = nullIndex
				sim.calls[ambulance.callIndex].numBumps += 1 # stats
				
				# consider new dispatch for bumped call
				addEvent!(sim.eventList; parentEvent = event, form = considerDispatch, time = sim.time, call = sim.calls[ambulance.callIndex])
			else
				ambWasOnRoute = false
			end
			
			# stats:
			if ambWasOnRoute
				ambulance.totalTravelDuration += sim.time - ambulance.route.startTime
				ambulance.totalTravelDistance += calcRouteDistance!(sim, ambulance.route, sim.time)
				if ambulance.status == ambGoingToCall
					ambulance.totalBusyDuration += sim.time - ambulance.route.startTime
				end
			end
			
			# dispatch chosen ambulance
			addEvent!(sim.eventList; parentEvent = event, form = ambDispatched, time = sim.time, ambulance = ambulance, call = call)
		end
		
################
	
	elseif eventForm == ambDispatched
		ambulance = sim.ambulances[event.ambIndex]
		status = ambulance.status # shorthand
		@assert(isFree(status) || status == ambGoingToCall)
		call = sim.calls[event.callIndex]
		@assert(call.status == callScreening || call.status == callQueued || call.status == callWaitingForAmb) # callWaitingForAmb if call bumped
		
		# stats:
		if status == ambIdleAtStation
			ambulance.numDispatchesFromStation += 1
		elseif isGoingToStation(status)
			ambulance.numDispatchesOnRoad += 1
		elseif status == ambGoingToCall
			ambulance.numDispatchesOnRoad += 1
			ambulance.numRedispatches += 1
		elseif status == ambFreeAfterCall
			ambulance.numDispatchesOnFree += 1
		else error()
		end
		
		setAmbStatus!(ambulance, ambGoingToCall, sim.time)
		# ambulance.stationIndex
		ambulance.callIndex = call.index
		changeRoute!(sim, ambulance.route, sim.responseTravelPriorities[call.priority], sim.time, call.location, call.nearestNodeIndex)
		
		call.status = callWaitingForAmb
		call.ambIndex = event.ambIndex
		call.dispatchTime = sim.time # stats
		call.queuedDuration = call.wasQueued ? sim.time - (call.arrivalTime + call.dispatchDelay) : 0.0
		call.ambDispatchLoc = ambulance.route.startLoc # same result as getRouteCurrentLocation!(sim.net, ambulance.route, sim.time)
		call.ambStatusBeforeDispatch = ambulance.prevStatus
		
		addEvent!(sim.eventList; parentEvent = event, form = ambReachesCall, time = ambulance.route.endTime, ambulance = ambulance, call = call)
		
		if sim.moveUpData.useMoveUp
			m = sim.moveUpData.moveUpModule # shorthand
			if m == compTableModule || m == ddsmModule || m == zhangIpModule || m == temp0Module || m == temp1Module || m == temp2Module
				addEvent!(sim.eventList; parentEvent = event, form = considerMoveUp, time = sim.time, ambulance = ambulance, addEventToAmb = false)
			end
		end
		
################
	
	elseif eventForm == ambReachesCall
		ambulance = sim.ambulances[event.ambIndex]
		@assert(ambulance.status == ambGoingToCall)
		call = sim.calls[event.callIndex]
		@assert(call.status == callWaitingForAmb)
		
		setAmbStatus!(ambulance, ambAtCall, sim.time)
		# ambulance.stationIndex
		# ambulance.callIndex
		ambulance.totalTravelDuration += sim.time - ambulance.route.startTime # stats
		ambulance.totalTravelDistance += calcRouteDistance!(sim, ambulance.route, sim.time) # stats
		ambulance.totalBusyDuration += sim.time - ambulance.route.startTime # stats
		ambulance.numCallsTreated += 1 # stats
		
		call.status = callOnSceneTreatment
		call.ambArrivalTime = sim.time # stats
		call.responseDuration = sim.time - call.arrivalTime # stats
		
		# transport call to hospital if needed, otherwise amb becomes free
		if call.transport
			# transport to hospital
			addEvent!(sim.eventList; parentEvent = event, form = ambGoesToHospital, time = sim.time + call.onSceneDuration, ambulance = ambulance, call = call)
		else
			# amb becomes free
			addEvent!(sim.eventList; parentEvent = event, form = ambBecomesFree, time = sim.time + call.onSceneDuration, ambulance = ambulance, call = call)
		end
		
################
	
	elseif eventForm == ambGoesToHospital
		ambulance = sim.ambulances[event.ambIndex]
		@assert(ambulance.status == ambAtCall)
		call = sim.calls[event.callIndex]
		@assert(call.status == callOnSceneTreatment)
		@assert(call.transport)
		
		ambulance.totalBusyDuration += call.onSceneDuration # stats
		
		# if hospital not specified for call, find closest hospital
		hospitalIndex = call.hospitalIndex
		if hospitalIndex == nullIndex
			hospitalIndex = nearestHospitalToCall!(sim, call, lowPriority)
		end
		@assert(hospitalIndex != nullIndex)
		hospital = sim.hospitals[hospitalIndex]
		
		setAmbStatus!(ambulance, ambGoingToHospital, sim.time)
		# ambulance.stationIndex
		# ambulance.callIndex
		changeRoute!(sim, ambulance.route, lowPriority, sim.time, hospital.location, hospital.nearestNodeIndex)
		
		call.status = callGoingToHospital
		call.hospitalIndex = hospitalIndex
		
		addEvent!(sim.eventList; parentEvent = event, form = ambReachesHospital, time = ambulance.route.endTime, ambulance = ambulance, call = call)
		
################
	
	elseif eventForm == ambReachesHospital
		ambulance = sim.ambulances[event.ambIndex]
		@assert(ambulance.status == ambGoingToHospital)
		call = sim.calls[event.callIndex]
		@assert(call.status == callGoingToHospital)
		@assert(call.transport)
		
		setAmbStatus!(ambulance, ambAtHospital, sim.time)
		# ambulance.stationIndex
		# ambulance.callIndex
		ambulance.totalTravelDuration += sim.time - ambulance.route.startTime # stats
		ambulance.totalTravelDistance += calcRouteDistance!(sim, ambulance.route, sim.time) # stats
		ambulance.totalBusyDuration += sim.time - ambulance.route.startTime # stats
		ambulance.numCallsTransported += 1 # stats
		
		call.status = callAtHospital
		call.hospitalArrivalTime = sim.time # stats
		sim.hospitals[call.hospitalIndex].numCalls += 1 # stats
		
		addEvent!(sim.eventList; parentEvent = event, form = ambBecomesFree, time = sim.time + call.handoverDuration, ambulance = ambulance, call = call)
		
################
	
	elseif eventForm == ambBecomesFree
		ambulance = sim.ambulances[event.ambIndex]
		@assert(ambulance.status == ambAtCall || ambulance.status == ambAtHospital)
		call = sim.calls[event.callIndex]
		@assert(call.status == callOnSceneTreatment || call.status == callAtHospital)
		
		# remove call, processing is finished
		delete!(sim.currentCalls, call)
		call.status = callProcessed
		
		setAmbStatus!(ambulance, ambFreeAfterCall, sim.time)
		ambulance.callIndex = nullIndex
		ambulance.totalBusyDuration += call.transport ? call.handoverDuration : call.onSceneDuration # stats
		# if call.transport == true, already added call.onSceneDuration to totalBusyDuration in ambGoesToHospital event
		
		# if queued call exists, respond
		# otherwise return to station
		if length(sim.queuedCallList) > 0
			call = getNextCall!(sim.queuedCallList)
			@assert(call != nothing)
			
			# dispatch ambulance
			addEvent!(sim.eventList; parentEvent = event, form = ambDispatched, time = sim.time, ambulance = ambulance, call = call)
		else
			# return to station (can be cancelled by move up)
			addEvent!(sim.eventList; parentEvent = event, form = ambReturnsToStation, time = sim.time, ambulance = ambulance)
			
			# consider move up
			# note that this event is created after above ambReturnsToStation event, so that it will happen first and can delete the ambReturnsToStation event if needed
			if sim.moveUpData.useMoveUp
				m = sim.moveUpData.moveUpModule # shorthand
				if m == compTableModule || m == dmexclpModule || m == priorityListModule || m == zhangIpModule || m == temp0Module || m == temp1Module || m == temp2Module
					addEvent!(sim.eventList; parentEvent = event, form = considerMoveUp, time = sim.time, ambulance = ambulance, addEventToAmb = false)
				end
			end
		end
		
################
	
	elseif eventForm == ambReturnsToStation
		ambulance = sim.ambulances[event.ambIndex]
		@assert(ambulance.status == ambFreeAfterCall)
		@assert(event.callIndex == nullIndex)
		
		station = sim.stations[ambulance.stationIndex]
		
		setAmbStatus!(ambulance, ambReturningToStation, sim.time)
		changeRoute!(sim, ambulance.route, lowPriority, sim.time, station.location, station.nearestNodeIndex)
		
		addEvent!(sim.eventList; parentEvent = event, form = ambReachesStation, time = ambulance.route.endTime, ambulance = ambulance)
		
################
	
	elseif eventForm == ambReachesStation
		ambulance = sim.ambulances[event.ambIndex]
		@assert(isGoingToStation(ambulance.status))
		@assert(event.callIndex == nullIndex)
		
		setAmbStatus!(ambulance, ambIdleAtStation, sim.time)
		# ambulance.stationIndex
		# ambulance.callIndex
		ambulance.totalTravelDuration += sim.time - ambulance.route.startTime # stats
		ambulance.totalTravelDistance += calcRouteDistance!(sim, ambulance.route, sim.time) # stats
		
		ambulance.event = Event() # no event currently
		
################
	
	elseif eventForm == considerMoveUp
		ambulance = sim.ambulances[event.ambIndex] # ambulance that triggered consideration of move up
		@assert(event.callIndex == nullIndex)
		mud = sim.moveUpData # shorthand
		@assert(mud.useMoveUp)
		
		# call move up function
		movableAmbs = [] # init
		ambStations = [] # init
		
		if sim.resim.use
			(movableAmbs, ambStations) = resimMoveUp(sim)
		else
			if mud.moveUpModule == compTableModule
				(movableAmbs, ambStations) = compTableMoveUp(sim)
			elseif mud.moveUpModule == ddsmModule
				(movableAmbs, ambStations) = ddsmMoveUp(sim)
			elseif mud.moveUpModule == dmexclpModule
				(movableAmbs, ambStations) = dmexclpMoveUp(sim, ambulance)
			elseif mud.moveUpModule == priorityListModule
				(movableAmbs, ambStations) = priorityListMoveUp(sim, ambulance)
			elseif mud.moveUpModule == zhangIpModule
				(movableAmbs, ambStations) = zhangIpMoveUp(sim)
			elseif mud.moveUpModule == temp0Module
				(movableAmbs, ambStations) = temp0MoveUp(sim)
			elseif mud.moveUpModule == temp1Module
				(movableAmbs, ambStations) = temp1MoveUp(sim)
			elseif mud.moveUpModule == temp2Module
				(movableAmbs, ambStations) = temp2MoveUp(sim)
			else
				error()
			end
		end
		
		for i = 1:length(movableAmbs)
			ambulance = movableAmbs[i]
			station = ambStations[i]
			
			# move up ambulance if ambulance station has changed
			if ambulance.stationIndex != station.index
				
				if isGoingToStation(ambulance.status) # ambulance.event.form == ambReachesStation
					# delete station arrival event for this ambulance
					deleteEvent!(sim.eventList, ambulance.event)
					
					ambulance.totalTravelDuration += sim.time - ambulance.route.startTime # stats
					ambulance.totalTravelDistance += calcRouteDistance!(sim, ambulance.route, sim.time) # stats
				elseif ambulance.status == ambFreeAfterCall && ambulance.event.form == ambReturnsToStation
					deleteEvent!(sim.eventList, ambulance.event)
				end
				
				ambulance.stationIndex = station.index
				
				addEvent!(sim.eventList; parentEvent = event, form = ambMoveUpToStation, time = sim.time, ambulance = ambulance)
			end
		end
		
################
	
	elseif eventForm == ambMoveUpToStation
		ambulance = sim.ambulances[event.ambIndex]
		@assert(event.callIndex == nullIndex)
		
		station = sim.stations[ambulance.stationIndex] # station to move up to
		
		# stats
		status = ambulance.status # shorthand
		if status == ambIdleAtStation
			ambulance.numMoveUpsFromStation += 1
		elseif isGoingToStation(status)
			ambulance.numMoveUpsOnRoad += 1
		elseif status == ambFreeAfterCall
			ambulance.numMoveUpsOnFree += 1
		else error()
		end
		
		setAmbStatus!(ambulance, ambMovingUpToStation, sim.time)
		ambulance.stationIndex = station.index
		ambulance.callIndex = nullIndex
		changeRoute!(sim, ambulance.route, lowPriority, sim.time, station.location, station.nearestNodeIndex)
		
		addEvent!(sim.eventList; parentEvent = event, form = ambReachesStation, time = ambulance.route.endTime, ambulance = ambulance)
		
################
	
	# elseif eventForm == ambRedirected
		# ambulance = sim.ambulances[event.ambIndex]
		# call = sim.calls[event.callIndex]
		# @assert(ambulance.status == ambGoingToCall)
		
		# # may alter ambDispatched event to deal with redirects
		
		# # not yet used
		# error()
		
################
	
	else
		# unspecified event
		error()
	end
	
end
