# load the resimulation data from the events file
# changes sim.resim.use to be false unless all checks are passed
function initResimulation!(sim::Simulation)
	resim = sim.resim # shorthand
	assert(resim.use)
	resim.use = false # until checks have passed
	
	eventsFilename = sim.outputFiles["events"].path
	println("resimulating, based on events from file: ", eventsFilename)
	
	if !isfile(eventsFilename)
		println("cannot resimulate, events file not found")
		return
	end
	
	# read events file
	(events, fileEnded, inputFiles, fileChecksums) = readEventsFile(eventsFilename)
	
	if !fileEnded
		println("cannot resimulate, events file closed before end")
		return
	end
	
	# check that checksum values of input files are same as in events file
	allMatch = true
	for i = 1:length(inputFiles)
		if fileChecksums[i] != sim.inputFiles[inputFiles[i]].checksum
			println(" checksum mismatch for file: ", inputFiles[i])
			allMatch = false
		end
	end
	if !allMatch
		println("cannot resimulate, input file checksums do not match those in events file")
		return
	end
	
	# all checks have passed, can resimulate
	resim.use = true
	resim.events = events
	resim.prevEventIndex = 0
	resim.timeTolerance = 1e-5 / 2 + 10*eps()
end

function resimCheckCurrentEvent!(sim::Simulation, event::Event)
	resim = sim.resim # shorthand
	assert(resim.use)
	
	resim.prevEventIndex += 1 # go to event after previous (should be current)
	
	resimEvent = resim.events[resim.prevEventIndex] # shorthand
	eventsMatch = true
	if abs(event.time - resimEvent.time) > resim.timeTolerance
		println("mismatching event time")
		eventsMatch = false
	elseif event.form != resimEvent.form
		println("mismatching event form")
		eventsMatch = false
	elseif event.ambIndex != resimEvent.ambIndex
		println("mismatching event ambulance index")
		eventsMatch = false
	elseif event.callIndex != resimEvent.callIndex
		println("mismatching event call index")
		eventsMatch = false
	elseif !(event.ambIndex == nullIndex && resimEvent.stationIndex == nullIndex || sim.ambulances[event.ambIndex].stationIndex == resimEvent.stationIndex)
		println("mismatching event station index")
		eventsMatch = false
	end
	if !eventsMatch
		@show event
		@show resimEvent
		error("resimulation event does not match current event")
	end
end

# select ambulance to dispatch to a call
# assumes ambulance mobilisation delay is 0
function resimNearestFreeAmbToCall(sim::Simulation, call::Call)
	resim = sim.resim # shorthand
	assert(resim.use)
	
	# check that previous event was to dispatch for this call
	resimEvent = resim.events[resim.prevEventIndex]
	assert(resimEvent.form == considerDispatch)
	assert(resimEvent.callIndex == call.index)
	
	resimEvent = resim.events[resim.prevEventIndex + 1] # go to next event, should be dispatch event (assuming ambulance mobilisation delay is 0)
	if resimEvent.form == ambDispatched
		assert(resimEvent.callIndex == call.index)
		assert(abs(sim.time - resimEvent.time) <= resim.timeTolerance)
		ambIndex = resimEvent.ambIndex
	else
		ambIndex = nullIndex
	end
	
	return ambIndex
end

# move up function for re-simulation
# assumes ambulance move up delay is 0
function resimMoveUp(sim::Simulation)
	resim = sim.resim # shorthand
	assert(resim.use)
	
	# check that previous event was to consider move up
	resimEvent = resim.events[resim.prevEventIndex]
	assert(resimEvent.form == considerMoveUp)
	
	# find all ambulances to move up
	movableAmbs = Vector{Ambulance}(0)
	ambStations = Vector{Station}(0)
	i = resim.prevEventIndex + 1 # next event, should be move up event (assuming ambulance move up delay is 0)
	resimEvent = resim.events[i]
	while abs(sim.time - resimEvent.time) <= resim.timeTolerance && resimEvent.form == ambMoveUp
		push!(movableAmbs, sim.ambulances[resimEvent.ambIndex])
		push!(ambStations, sim.stations[resimEvent.stationIndex])
		
		# go to next event, multiple ambulances may be involved in a single move up
		i += 1
		resimEvent = resim.events[i]
	end
	
	# reverse order of move up events...
	movableAmbs = flipdim(movableAmbs, 1)
	ambStations = flipdim(ambStations, 1)
	
	return movableAmbs, ambStations
end
