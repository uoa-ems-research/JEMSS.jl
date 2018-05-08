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
	i = findlast(e -> e.time >= event.time, eventList) + 1
	insert!(eventList, i, event)
	
	if checkMode
		# check time ordering of events
		assert(issorted(eventList, by = e -> e.time, rev = true))
	end
	
	return event
end

# for adding call arrival to event list
function addEvent!(eventList::Vector{Event}, call::Call)
	addEvent!(eventList; form = callArrives, time = call.arrivalTime, call = call)
end

# get next event from list
function getNextEvent!(eventList::Vector{Event})
	return length(eventList) > 0 ? pop!(eventList) : Event()
end

function printEvent(event::Event)
	println("Event:")
	println("form = ", string(event.form))
	println("time = ", event.time)
	println("ambIndex = ", event.ambIndex)
	println("callIndex = ", event.callIndex)
end

# delete event in eventList
function deleteEvent!(eventList::Vector{Event}, event::Event)
	i = findfirst(e -> e == event, eventList)
	assert(i != 0)
	assert(findnext(e -> e == event, eventList, i + 1) == 0)
	deleteat!(eventList, i)
end
