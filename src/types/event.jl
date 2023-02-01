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

# add new event to list, maintain sorting by time
function addEvent!(eventList::Vector{Event};
    parentEvent::Union{Event,Nothing}=nothing, form::EventForm=nullEvent, time::Float=nullTime,
    ambulance::Union{Ambulance,Nothing}=nothing, call::Union{Call,Nothing}=nothing, station::Union{Station,Nothing}=nothing,
    addEventToAmb::Bool=true)

    event = Event()
    event.parentIndex = parentEvent === nothing ? nullIndex : parentEvent.index
    event.form = form
    event.time = time
    event.ambIndex = ambulance === nothing ? nullIndex : ambulance.index
    event.callIndex = call === nothing ? nullIndex : call.index
    event.stationIndex = station === nothing ? nullIndex : station.index

    # set next event of ambulance
    if addEventToAmb && ambulance !== nothing
        ambulance.event = event
    end

    # find where to insert event into list
    # maintain sorting by time, events nearer to end of list are sooner
    i = something(findlast(e -> e.time >= event.time, eventList), 0) + 1
    insert!(eventList, i, event)

    if checkMode
        # check time ordering of events
        for i = 1:length(eventList)-1
            @assert(eventList[i].time >= eventList[i+1].time)
        end
        # @assert(issorted(eventList, by = e -> e.time, rev = true)) # slow
    end

    return event
end

# for adding call arrival to event list
function addEvent!(eventList::Vector{Event}, call::Call)
    addEvent!(eventList; form=callArrives, time=call.arrivalTime, call=call)
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
    i = something(findfirst(isequal(event), eventList))
    @assert(i !== nothing)
    @assert(findnext(e -> e == event, eventList, i + 1) === nothing)
    deleteat!(eventList, i)
end
