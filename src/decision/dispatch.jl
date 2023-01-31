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

# return true if ambulance can be dispatched to a call
function isAmbDispatchable(sim::Simulation, ambulance::Ambulance, call::Call)
    status = ambulance.status
    if isFree(status)
        return true
    elseif in(status, (ambMobilising, ambGoingToCall))
        return isAmbRedispatchable(sim, ambulance, sim.calls[ambulance.callIndex], call)
    end
    return false
end

# return true if ambulance may be redispatched from its current call to a different call, false otherwise.
function isAmbRedispatchable(sim::Simulation, ambulance::Ambulance, fromCall::Call, toCall::Call)
    @assert(in(ambulance.status, (ambMobilising, ambGoingToCall)))
    @assert(ambulance.callIndex == fromCall.index)
    sim.redispatch.allow || return false
    return sim.redispatch.conditions[Int(fromCall.priority), Int(toCall.priority)]
end

# return the index of the nearest dispatchable ambulance for a call
function findNearestDispatchableAmb!(sim::Simulation, call::Call)

    # nearest node to call, this is independent of chosen ambulance
    (node2, dist2) = (call.nearestNodeIndex, call.nearestNodeDist)

    # for ambulances that can be dispatched, find the one with shortest travel time to call
    ambIndex = nullIndex # nearest free ambulance
    minDuration = Inf
    for amb in sim.ambulances
        if isAmbDispatchable(sim, amb, call)
            # get time that ambulance will mobilise
            mobilisationTime = sim.time
            if sim.mobilisationDelay.use
                if amb.status == ambIdleAtStation
                    mobilisationTime += sim.mobilisationDelay.expectedDuration
                elseif amb.status == ambMobilising
                    # dispatcher does not have perfect knowledge of actual dispatch duration
                    mobilisationTime += max(0, sim.mobilisationDelay.expectedDuration - (sim.time - amb.dispatchTime))
                end
            end
            responseDuration = mobilisationTime - sim.time # will add travel time next

            # duration on network
            travelMode = getTravelMode!(sim.travel, sim.responseTravelPriorities[call.priority], sim.time; startTime=mobilisationTime)
            (node1, time1) = getRouteNextNode!(sim, amb.route, travelMode.index, sim.time) # next/nearest node in ambulance route
            responseDuration += shortestPathTravelTime(sim.net, travelMode.index, node1, node2) # time spent on network
            time2 = offRoadTravelTime(travelMode, dist2) # time to reach nearest node
            responseDuration += time1 + time2 # add time to get on and off network

            if minDuration > responseDuration
                ambIndex = amb.index
                minDuration = responseDuration
            end
        end
    end

    return ambIndex
end
findNearestFreeAmbToCall!(sim::Simulation, call::Call) = findNearestDispatchableAmb!(sim, call) # compat
