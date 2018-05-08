# run simulation, show progress whenever system time increases by timeStep (seconds)
function simulate!(sim::Simulation; timeStep::Float = 1.0)
	println("running simulation...")
	startTime = time()
	nextTime = startTime + timeStep
	printProgress() = print(@sprintf("\rsim duration: %-9.2f real duration: %.2f seconds", sim.time - sim.startTime, time()-startTime))
	while !sim.complete
		simulateNextEvent!(sim)
		if time() > nextTime
			printProgress()
			nextTime += timeStep
		end
	end
	printProgress()
	println("\n...simulation complete")
end

# simulate up until given time
function simulateToTime!(sim::Simulation, time::Float)
	while !sim.complete && sim.eventList[end].time <= time # eventList must remain sorted by non-increasing time
		simulateNextEvent!(sim)
	end
end

function simulateToEnd!(sim::Simulation)
	simulateToTime!(sim, Inf)
end

# set sim.backup to copy of sim (before running)
# note that for reducing memory usage, sim.backup does not contain: net, travel, grid, resim
function backupSim!(sim::Simulation)
	assert(!sim.used)
	
	# remove net, travel, grid, and resim from sim before copying sim
	(net, travel, grid, resim) = (sim.net, sim.travel, sim.grid, sim.resim)
	(sim.net, sim.travel, sim.grid, sim.resim) = (Network(), Travel(), Grid(), Resimulation())
	
	sim.backup = deepcopy(sim)
	
	(sim.net, sim.travel, sim.grid, sim.resim) = (net, travel, grid, resim)
end

# reset sim from sim.backup
function resetSim!(sim::Simulation)
	assert(!sim.backup.used)
	
	if sim.used
		resetCalls!(sim) # reset calls from sim.backup, need to do this before resetting sim.time
		
		fnames = Set(fieldnames(sim))
		fnamesDontCopy = Set([:backup, :net, :travel, :grid, :resim, :calls]) # will not (yet) copy these fields from sim.backup to sim
		# note that sim.backup does not contain net, travel, grid, or resim
		setdiff!(fnames, fnamesDontCopy) # remove fnamesDontCopy from fnames
		for fname in fnames
			try
				setfield!(sim, fname, deepcopy(getfield(sim.backup, fname)))
			end
		end
		
		# reset resimulation state
		sim.resim.prevEventIndex = 0
		
		# reset travel state
		sim.travel.recentSetsStartTimesIndex = 1
	end
end

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
		assert(sim.endTime == nullTime)
		assert(sim.complete == false)
		sim.endTime = sim.time
		sim.complete = true
	end
end

# simulate event
function simulateEvent!(sim::Simulation, event::Event)
	# format:
	# next event may change relevant ambulance / call fields at event.time
	# event may then trigger future events / cancel scheduled events
	
	assert(sim.time == event.time)
	
	eventForm = event.form
	if eventForm == nullEvent
		error("null event")
		
################

	elseif eventForm == ambGoesToSleep
		ambulance = sim.ambulances[event.ambIndex]
		assert(ambulance.status == ambIdleAtStation) # ambulance should be at station before sleeping
		assert(event.callIndex == nullIndex)
	
		ambulance.status = ambSleeping
		# ambulance.stationIndex
		# ambulance.callIndex
		
		addEvent!(sim.eventList; parentEvent = event, form = ambWakesUp, time = sim.time + sleepDuration, ambulance = ambulance)
		
################

	elseif eventForm == ambWakesUp
		ambulance = sim.ambulances[event.ambIndex]
		assert(ambulance.status == ambSleeping) # ambulance should have been sleeping
		assert(event.callIndex == nullIndex)
		
		ambulance.status = ambIdleAtStation
		# ambulance.stationIndex
		# ambulance.callIndex
		
################

	elseif eventForm == callArrives
		assert(event.ambIndex == nullIndex) # no ambulance should be assigned yet
		call = sim.calls[event.callIndex]
		assert(call.status == callNullStatus)
		
		push!(sim.currentCalls, sim.calls[event.callIndex])
		call.status = callScreening
		
		addEvent!(sim.eventList; parentEvent = event, form = considerDispatch, time = sim.time + call.dispatchDelay, call = call)
		
		# add next call arrival to event queue
		if call.index < length(sim.calls)
			addEvent!(sim.eventList, sim.calls[call.index + 1])
		end
		
################

	elseif eventForm == considerDispatch
		assert(event.ambIndex == nullIndex) # no ambulance should be assigned yet
		call = sim.calls[event.callIndex]
		assert(call.status == callScreening || call.status == callWaitingForAmb) # callWaitingForAmb if call bumped
		
		# find ambulance to respond
		if sim.resim.use
			ambIndex = resimFindAmbToDispatch(sim, call)
		else
			ambIndex = sim.findAmbToDispatch!(sim, call)
		end
		# if no ambulance found, add call to queue,
		# will respond when an ambulance becomes idle
		if ambIndex == nullIndex
			sim.addCallToQueue!(sim.queuedCallList, call)
			call.status = callQueued
			call.wasQueued = true # stats
		else
			ambulance = sim.ambulances[ambIndex]
			
			# remove any previously scheduled events
			ambWasOnRoute = true
			if ambulance.status == ambGoingToStation
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
				ambulance.totalTravelTime += sim.time - ambulance.route.startTime
				if ambulance.status == ambGoingToCall
					ambulance.totalBusyTime += sim.time - ambulance.route.startTime 
					ambulance.numDiversions += 1
				end
			end
			
			# dispatch chosen ambulance
			addEvent!(sim.eventList; parentEvent = event, form = ambDispatched, time = sim.time, ambulance = ambulance, call = call)
		end
		
################

	elseif eventForm == ambDispatched
		ambulance = sim.ambulances[event.ambIndex]
		status = ambulance.status # shorthand
		assert(status == ambIdleAtStation || status == ambGoingToCall ||
			status == ambGoingToStation || status == ambAtCall || status == ambAtHospital)
		call = sim.calls[event.callIndex]
		assert(call.status == callScreening || call.status == callQueued || call.status == callWaitingForAmb) # callWaitingForAmb if call bumped
		
		# stats:
		if status == ambIdleAtStation
			ambulance.atStationDispatches += 1
		elseif status == ambGoingToCall || status == ambGoingToStation
			ambulance.onRoadDispatches += 1
		elseif status == ambAtCall || status == ambAtHospital
			ambulance.afterServiceDispatches += 1
		end
		
		ambulance.status = ambGoingToCall
		# ambulance.stationIndex
		ambulance.callIndex = call.index
		changeRoute!(sim, ambulance.route, call.priority, sim.time, call.location, call.nearestNodeIndex)
		
		call.status = callWaitingForAmb
		call.ambIndex = event.ambIndex
		call.dispatchTime = sim.time # stats
		
		addEvent!(sim.eventList; parentEvent = event, form = ambReachesCall, time = ambulance.route.endTime, ambulance = ambulance, call = call)
		
		if sim.moveUpData.useMoveUp
			m = sim.moveUpData.moveUpModule # shorthand
			if m == compTableModule || m == zhangIpModule || m == temp1Module || m == temp2Module
				addEvent!(sim.eventList; parentEvent = event, form = considerMoveUp, time = sim.time, ambulance = ambulance, addEventToAmb = false)
			end
		end
		
################

	elseif eventForm == ambReachesCall
		ambulance = sim.ambulances[event.ambIndex]
		assert(ambulance.status == ambGoingToCall)
		call = sim.calls[event.callIndex]
		assert(call.status == callWaitingForAmb)
		
		ambulance.status = ambAtCall
		# ambulance.stationIndex
		# ambulance.callIndex
		ambulance.totalTravelTime += sim.time - ambulance.route.startTime # stats
		ambulance.totalBusyTime += sim.time - ambulance.route.startTime # stats
		ambulance.numCallsTreated += 1 # stats
		
		call.status = callOnSceneCare
		call.ambArrivalTime = sim.time # stats
		call.responseTime = sim.time - call.arrivalTime # stats
		
		# transport call to hospital if needed, otherwise amb becomes idle
		if call.transfer
			# transport to hospital
			addEvent!(sim.eventList; parentEvent = event, form = ambGoesToHospital, time = sim.time + call.onSceneDuration, ambulance = ambulance, call = call)
		else
			# amb becomes idle
			addEvent!(sim.eventList; parentEvent = event, form = ambBecomesIdle, time = sim.time + call.onSceneDuration, ambulance = ambulance, call = call)
		end
		
################

	elseif eventForm == ambGoesToHospital
		ambulance = sim.ambulances[event.ambIndex]
		assert(ambulance.status == ambAtCall)
		call = sim.calls[event.callIndex]
		assert(call.status == callOnSceneCare)
		
		# if hospital not specified for call, find closest hospital
		hospitalIndex = call.hospitalIndex
		if hospitalIndex == nullIndex
			travelMode = getTravelMode!(sim.travel, lowPriority, sim.time)
			hospitalIndex = nearestHospitalToCall(travelMode, call)
		end
		assert(hospitalIndex != nullIndex)
		hospital = sim.hospitals[hospitalIndex]
		
		ambulance.status = ambGoingToHospital
		# ambulance.stationIndex
		# ambulance.callIndex
		changeRoute!(sim, ambulance.route, lowPriority, sim.time, hospital.location, hospital.nearestNodeIndex)
		
		call.status = callGoingToHospital
		call.hospitalIndex = hospitalIndex
		
		addEvent!(sim.eventList; parentEvent = event, form = ambReachesHospital, time = ambulance.route.endTime, ambulance = ambulance, call = call)
		
################

	elseif eventForm == ambReachesHospital
		ambulance = sim.ambulances[event.ambIndex]
		assert(ambulance.status == ambGoingToHospital)
		call = sim.calls[event.callIndex]
		assert(call.status == callGoingToHospital)

		ambulance.status = ambAtHospital
		# ambulance.stationIndex
		# ambulance.callIndex
		ambulance.totalTravelTime += sim.time - ambulance.route.startTime # stats
		ambulance.totalBusyTime += sim.time - ambulance.route.startTime # stats
		ambulance.numCallsTransferred += 1 # stats
		
		call.status = callAtHospital
		call.hospitalArrivalTime = sim.time # stats
		sim.hospitals[call.hospitalIndex].numTransfers += 1 # stats
		
		addEvent!(sim.eventList; parentEvent = event, form = ambBecomesIdle, time = sim.time + call.transferDuration, ambulance = ambulance, call = call)
		
################

	elseif eventForm == ambBecomesIdle
		ambulance = sim.ambulances[event.ambIndex]
		assert(ambulance.status == ambAtCall || ambulance.status == ambAtHospital)
		call = sim.calls[event.callIndex]
		assert(call.status == callOnSceneCare || call.status == callAtHospital)
		
		# remove call, processing is finished
		delete!(sim.currentCalls, call)
		call.status = callProcessed
		
		ambulance.totalBusyTime += call.onSceneDuration + call.transfer * call.transferDuration # stats
		
		# if queued call exists, respond
		# otherwise return to station
		if length(sim.queuedCallList) > 0
			call = getNextCall!(sim.queuedCallList)
			assert(call != nothing)

			# dispatch ambulance
			addEvent!(sim.eventList; parentEvent = event, form = ambDispatched, time = sim.time, ambulance = ambulance, call = call)
		else
			station = sim.stations[ambulance.stationIndex]
			
			ambulance.status = ambGoingToStation
			# ambulance.stationIndex
			ambulance.callIndex = nullIndex
			changeRoute!(sim, ambulance.route, lowPriority, sim.time, station.location, station.nearestNodeIndex)
			
			addEvent!(sim.eventList; parentEvent = event, form = ambReachesStation, time = ambulance.route.endTime, ambulance = ambulance)
			
			if sim.moveUpData.useMoveUp
				m = sim.moveUpData.moveUpModule # shorthand
				if m == compTableModule || m == dmexclpModule || m == priorityListModule || m == zhangIpModule || m == temp1Module || m == temp2Module
					addEvent!(sim.eventList; parentEvent = event, form = considerMoveUp, time = sim.time, ambulance = ambulance, addEventToAmb = false)
				end
			end
		end
		
################

	elseif eventForm == ambReachesStation
		ambulance = sim.ambulances[event.ambIndex]
		assert(ambulance.status == ambGoingToStation)
		assert(event.callIndex == nullIndex)
		
		ambulance.status = ambIdleAtStation
		# ambulance.stationIndex
		# ambulance.callIndex
		ambulance.totalTravelTime += sim.time - ambulance.route.startTime # stats
		
		ambulance.event = Event() # no event currently
		
################

	elseif eventForm == considerMoveUp
		ambulance = sim.ambulances[event.ambIndex] # ambulance that triggered consideration of move up
		assert(event.callIndex == nullIndex)
		mud = sim.moveUpData # shorthand
		assert(mud.useMoveUp)
		
		# call move up function
		movableAmbs = [] # init
		ambStations = [] # init
		
		if sim.resim.use
			(movableAmbs, ambStations) = resimMoveUp(sim)
		else
			if mud.moveUpModule == compTableModule
				(movableAmbs, ambStations) = compTableMoveUp(sim)
			elseif mud.moveUpModule == dmexclpModule
				(movableAmbs, ambStations) = dmexclpMoveUp(sim, ambulance)
			elseif mud.moveUpModule == priorityListModule
				(movableAmbs, ambStations) = priorityListMoveUp(sim, ambulance)
			elseif mud.moveUpModule == zhangIpModule
				(movableAmbs, ambStations) = zhangIpMoveUp(sim)
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
				
				if ambulance.status == ambGoingToStation # ambulance.event.form == ambReachesStation
					# delete station arrival event for this ambulance
					deleteEvent!(sim.eventList, ambulance.event)
					
					ambulance.totalTravelTime += sim.time - ambulance.route.startTime # stats
				end
				
				# ambulance.status
				ambulance.stationIndex = station.index
				# ambulance.callIndex
				
				addEvent!(sim.eventList; parentEvent = event, form = ambMoveUp, time = sim.time, ambulance = ambulance)
			end
		end
		
################

	elseif eventForm == ambMoveUp
		ambulance = sim.ambulances[event.ambIndex]
		assert(event.callIndex == nullIndex)
		
		station = sim.stations[ambulance.stationIndex] # station to move up to
		
		ambulance.status = ambGoingToStation
		ambulance.stationIndex = station.index
		ambulance.callIndex = nullIndex
		changeRoute!(sim, ambulance.route, lowPriority, sim.time, station.location, station.nearestNodeIndex)
		
		addEvent!(sim.eventList; parentEvent = event, form = ambReachesStation, time = ambulance.route.endTime, ambulance = ambulance)
		
################

	# elseif eventForm == ambRedirected
		# ambulance = sim.ambulances[event.ambIndex]
		# call = sim.calls[event.callIndex]
		# assert(ambulance.status == ambGoingToCall)
		
		# # may alter ambDispatched event to deal with redirects
		
		# # not yet used
		# error()
		
################

	else
		# unspecified event
		error()
	end
	
end

# find nearest hospital to call location, given the travel mode to use
function nearestHospitalToCall(travelMode::TravelMode, call::Call)
	return travelMode.fNetTravel.fNodeNearestHospitalIndex[call.nearestNodeIndex]
end
