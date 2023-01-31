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

# get next call from list
function getNextCall!(queuedCallList::Vector{Call})
    return length(queuedCallList) > 0 ? pop!(queuedCallList) : nothing
end

# reset calls
function reset!(calls::Vector{Call})
    # do not reset fixed call parameters
    fnamesFixed = (:index, :priority, :transport, :hospitalIndex, :location,
        :arrivalTime, :dispatchDelay, :onSceneDuration, :handoverDuration,
        :nearestNodeIndex, :nearestNodeDist)
    fnames = setdiff(fieldnames(Call), fnamesFixed)

    # reset calls
    nullCall = Call()
    recentCallIndex = something(findlast(call -> call.status != callNullStatus, calls), 0) # index of last call that needs reset
    for fname in fnames
        val = getfield(nullCall, fname)
        if isprimitivetype(fieldtype(Call, fname))
            for i = 1:recentCallIndex
                setfield!(calls[i], fname, val)
            end
        elseif fieldtype(Call, fname) == Location
            for i = 1:recentCallIndex
                copy!(getfield(calls[i], fname), val) # faster than deepcopy for Location
            end
        else
            for i = 1:recentCallIndex
                setfield!(calls[i], fname, deepcopy(val))
            end
        end
    end
end
resetCalls!(sim::Simulation) = reset!(sim.calls) # compat

function setCallStatus!(call::Call, status::CallStatus, time::Float)
    @assert(call.statusSetTime <= time)

    # stats for previous status
    statusDuration = time - call.statusSetTime
    prevStatus = call.status # shorthand
    if prevStatus == callQueued
        call.queuedDuration += statusDuration
    elseif prevStatus == callWaitingForAmb
        call.waitingForAmbDuration += statusDuration
        call.ambGoingToCallDuration = statusDuration # overwrite previous value if have attempted dispatch to call multiple times
    elseif prevStatus == callGoingToHospital
        @assert(call.transportDuration == 0.0) # should not have been calculated yet
        call.transportDuration = statusDuration
    end

    # stats for new status
    if status == callQueued
        call.wasQueued = true
    elseif status == callWaitingForAmb
        if prevStatus == callScreening
            @assert(isapprox(time, call.arrivalTime + call.dispatchDelay)) # if dispatching after call screening, only delay should be dispatch delay
        end
        call.dispatchTime = time
    elseif status == callOnSceneTreatment
        @assert(call.ambArrivalTime == nullTime) # value should not have been set yet
        call.ambArrivalTime = time
        call.responseDuration = time - call.arrivalTime
        call.bumpedDuration = call.waitingForAmbDuration - call.ambGoingToCallDuration
    elseif status == callAtHospital
        call.hospitalArrivalTime = time
    elseif status == callProcessed
        @assert(call.serviceDuration == 0.0) # value should not have been set yet
        call.serviceDuration = time - call.arrivalTime
    end

    call.status = status
    call.statusSetTime = time
end

function initCalls!(sim::Simulation, calls::Vector{Call})
    for (i, call) in enumerate(calls)
        call.index = i

        # find nearest nodes
        if call.nearestNodeIndex == nullIndex
            (call.nearestNodeIndex, call.nearestNodeDist) = findNearestNode(sim.map, sim.grid, sim.net.fGraph.nodes, call.location)
        end
    end
end

# Set the calls of the simulation.
# Note that the old sim.calls will not be kept.
# I tried to write this function to allow sim.startTime to be changed,
# but that would also require sim stats capturing to change, so startTime cannot be changed.
# Mutates: sim, calls
function setSimCalls!(sim::Simulation, calls::Vector{Call})
    @assert(length(calls) >= 1)
    @assert(sim.startTime <= calls[1].arrivalTime)

    reset!(sim) # for safety, as changing calls part-way through a simulation could be tricky

    sim.numCalls = length(calls)

    initCalls!(sim, calls)

    # set first event to match first call arrival time; hacky
    firstCallEvents = filter(e -> e.form == callArrives, sim.eventList)
    @assert(length(firstCallEvents) == 1)
    firstCallEvents[1].time = calls[1].arrivalTime

    if isdefined(sim, :backup)
        sim.backup.numCalls = sim.numCalls
        sim.backup.eventList = deepcopy(sim.eventList)
    end

    reset!(calls)

    sim.calls = calls
end
