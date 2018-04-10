# add new event to list, maintain sorting by time
function addEvent!(eventList::Vector{Event};
	parentEvent::Event = Event(), form::EventForm = nullEvent, time::Float = nullTime, ambulance::Ambulance = Ambulance(), call::Call = Call(), addEventToAmb::Bool = true)
	
	event = Event()
	event.parentIndex = parentEvent.index
	event.form = form
	event.time = time
	event.ambIndex = ambulance.index
	event.callIndex = call.index
	
	# set next event of ambulance
	if addEventToAmb
		ambulance.event = event
	end
	
	# find where to insert event into list
	# maintain sorting by time, events nearer to end of list are sooner
	i = length(eventList) + 1
	while i > 1 && eventList[i-1].time < event.time
		i -= 1
	end
	insert!(eventList, i, event)
	
	if checkMode
		# check time ordering of events
		for i = 1:length(eventList)-1
			assert(eventList[i].time >= eventList[i+1].time)
		end
	end
	
	return event
end

# for adding call arrival to event list
function addEvent!(eventList::Vector{Event}, call::Call)
	addEvent!(eventList; form = callArrives, time = call.arrivalTime, call = call)
end

# get next event from list
function getNextEvent!(eventList::Vector{Event})
	nextEvent = Event()
	if length(eventList) > 0
		nextEvent = pop!(eventList)
	end
	return nextEvent
end

function printEvent(event::Event)
	println("Event:")
	println("form = ", string(event.form))
	println("time = ", event.time)
	println("ambIndex = ", event.ambIndex)
	println("callIndex = ", event.callIndex)
end

# delete event in eventList
# return true if matching event deleted from eventList, false otherwise
function deleteEvent!(eventList::Vector{Event}, event::Event)
	n = length(eventList)
	
	# find index of event in list
	index = 0
	numMatches = 0
	for i = 1:n
		if event == eventList[i]
			index = i
			numMatches += 1
		end
	end
	if numMatches == 0
		return false
	end
	assert(numMatches == 1) # should only have one matching event
	
	# remove matching event from eventList
	deleteat!(eventList, index)
	
	return true
end
