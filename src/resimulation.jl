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

# load the resimulation data from the events file
# changes sim.resim.use to be false unless all checks are passed
function initResim!(sim::Simulation; doPrint::Bool = false)
	resim = sim.resim # shorthand
	@assert(resim.use)
	resim.use = false # until checks have passed
	
	eventsFilename = sim.outputFiles["events"].path
	doPrint && println("resimulating, based on events from file: ", eventsFilename)
	
	if !isfile(eventsFilename)
		doPrint && println("cannot resimulate, events file not found")
		return
	end
	
	# read events file
	(events, eventsChildren, eventFilter, fileEnded, inputFiles, fileChecksums) = readEventsFile(eventsFilename)
	
	if !fileEnded
		doPrint && println("cannot resimulate, events file closed before end")
		return
	end
	
	# check that checksum values of input files are same as in events file
	allMatch = true
	for i = 1:length(inputFiles)
		if inputFiles[i] == "statsControl" continue end # allow this file to change, should not affect sim events
		if fileChecksums[i] != sim.inputFiles[inputFiles[i]].checksum
			doPrint && println(" checksum mismatch for file: ", inputFiles[i])
			allMatch = false
		end
	end
	if !allMatch
		doPrint && println("cannot resimulate, input file checksums do not match those in events file")
		return
	end
	
	# all checks have passed, can resimulate
	resim.use = true
	resim.events = events
	resim.eventsChildren = eventsChildren
	resim.eventFilter = eventFilter
	resim.doDispatch = eventFilter[ambDispatched]
	resim.doMoveUp = eventFilter[ambMoveUpToStation]
	resim.prevEventIndex = 0
	resim.timeTolerance = 1e-5 / 2
end

function resimCheckCurrentEvent!(sim::Simulation, event::Event)
	resim = sim.resim # shorthand
	@assert(resim.use)
	
	i = resim.prevEventIndex # shorthand
	if i == length(resim.events) || resim.events[i+1].index > event.index
		return # have run out of events to resimulate, or do not have current event to resimulate with
	end
	resim.prevEventIndex += 1 # go to event after previous (should be current)
	
	resimEvent = resim.events[resim.prevEventIndex] # shorthand
	eventsMatch = true
	@assert(resimEvent.index == event.index)
	if !isapprox(event.time, resimEvent.time; rtol = eps(Float), atol = resim.timeTolerance)
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
	elseif event.stationIndex != resimEvent.stationIndex
		println("mismatching event station index")
		eventsMatch = false
	end
	if !eventsMatch
		@show event
		@show resimEvent
		error("resimulation event does not match current event")
	end
end

# find ambulance to dispatch to a call
# assumes ambulance mobilisation delay is 0
function resimFindAmbToDispatch(sim::Simulation, call::Call)
	resim = sim.resim # shorthand
	@assert(resim.use && resim.doDispatch)
	eventIndex = sim.eventIndex # shorthand
	
	resimEvent = resim.events[resim.prevEventIndex]
	if resimEvent.index == eventIndex
		# check that previous event was to dispatch for this call
		@assert(resimEvent.form == considerDispatch)
		@assert(resimEvent.callIndex == call.index)
	end
	
	# find child event for dispatch
	ambIndex = nullIndex
	if length(resim.eventsChildren) < eventIndex return ambIndex end
	eventChildren = resim.eventsChildren[eventIndex]
	dispatchEvents = filter(event -> event.form == ambDispatched, eventChildren)
	@assert(length(dispatchEvents) <= 1)
	if length(dispatchEvents) == 1
		ambIndex = dispatchEvents[1].ambIndex
	end
	
	return ambIndex
end

# move up function for re-simulation
# assumes ambulance move up delay is 0
function resimMoveUp(sim::Simulation)
	resim = sim.resim # shorthand
	@assert(resim.use && resim.doMoveUp)
	eventIndex = sim.eventIndex # shorthand
	
	resimEvent = resim.events[resim.prevEventIndex]
	if resimEvent.index == eventIndex
		# check that previous event was to consider move up
		@assert(resimEvent.form == considerMoveUp)
	end
	
	# find all ambulances to move up
	movableAmbs = Vector{Ambulance}()
	ambStations = Vector{Station}()
	if length(resim.eventsChildren) < eventIndex return movableAmbs, ambStations end
	eventChildren = resim.eventsChildren[eventIndex]
	for event in reverse(eventChildren)
		@assert(event.form == ambMoveUpToStation)
		push!(movableAmbs, sim.ambulances[event.ambIndex])
		push!(ambStations, sim.stations[event.stationIndex])
	end
	
	return movableAmbs, ambStations
end
