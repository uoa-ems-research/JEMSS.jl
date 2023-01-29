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

# priority list, for newly freed ambulances only
# (a priority list which uses all free ambulances can be described as a nested compliance table)
# From thesis: "Simulation optimisation and Markov models for dynamic ambulance redeployment"

# initialise data relevant to move up
function initPriorityList!(sim::Simulation, priorityListFilename::String)
    # read priority list
    priorityList = readPriorityListFile(priorityListFilename)
    initPriorityList!(sim, priorityList)
end

# initialise data relevant to move up
function initPriorityList!(sim::Simulation, priorityList::PriorityList)
    # shorthand names:
    @unpack ambulances, numStations = sim
    pld = sim.moveUpData.priorityListData

    checkPriorityList(priorityList, sim)
    pld.priorityList = priorityList

    pld.stationNumFreeAmbs = Vector{Int}(undef, numStations)

    # check that ambulance to station assignments follow priority list
    pld.stationNumFreeAmbs[:] .= 0
    for i = 1:sim.numAmbs
        pld.stationNumFreeAmbs[ambulances[i].stationIndex] += 1
    end
    if !all(pld.stationNumFreeAmbs[i] == sum(pld.priorityList .== i) for i = 1:numStations)
        @warn("Number of ambulances at each station does not match priority list")
    end
end

function priorityListMoveUp(sim::Simulation, newlyFreedAmb::Ambulance)
    @assert(sim.moveUpData.useMoveUp)

    # shorthand:
    pld = sim.moveUpData.priorityListData
    @unpack priorityList, stationNumFreeAmbs = pld
    @unpack ambulances, stations, numAmbs = sim

    # calculate the number of free ambulances at (or travelling to) each station
    stationNumFreeAmbs[:] .= 0
    for i = 1:numAmbs
        # do not count newly freed ambulance, it has not been assigned a station
        if isAmbMovable(ambulances[i]) && i != newlyFreedAmb.index
            stationNumFreeAmbs[ambulances[i].stationIndex] += 1
        end
    end
    numFreeAmbs = sum(stationNumFreeAmbs) + 1

    # return first item in priority list that is not already filled
    stationIndex = nullIndex
    for i = 1:numFreeAmbs
        j = priorityList[i] # station index
        stationNumFreeAmbs[j] -= 1
        if stationNumFreeAmbs[j] < 0
            stationIndex = j
            break
        end
    end
    @assert(stationIndex != nullIndex)

    return [newlyFreedAmb], [stations[stationIndex]]
end

# check that priority list is valid
function checkPriorityList(priorityList::PriorityList;
    numAmbs::Int=nullIndex, numStations::Int=nullIndex,
    stationCapacities::Union{Vector{Int},Nothing}=nothing)

    m = length(priorityList) # number of ambulances
    @assert(numAmbs == nullIndex || numAmbs == m)
    @assert(all(v -> v >= 1, priorityList)) # need positive integers (station indices)
    @assert(numStations == nullIndex || maximum(priorityList) <= numStations)
    if stationCapacities != nothing
        n = length(stationCapacities) # number of stations
        @assert(numStations == nullIndex || numStations == n)
        # check that priority list does not violate station capacity constraints
        for j = 1:n
            @assert(count(isequal(j), priorityList) <= stationCapacities[j])
        end
    end
    return true
end
function checkPriorityList(priorityList::PriorityList, sim::Simulation)
    checkPriorityList(priorityList; numAmbs=sim.numAmbs, numStations=sim.numStations, stationCapacities=[s.capacity for s in sim.stations])
end

# returns a randomly generated priority list
function makeRandPriorityList(numAmbs::Int, numStations::Int;
    stationCapacities::Union{Vector{Int},Nothing}=nothing,
    rng::AbstractRNG=GLOBAL_RNG)::PriorityList

    @assert(numStations > 0)
    if stationCapacities == nothing
        return rand(rng, 1:numStations, numAmbs)
    else
        @assert(numStations == length(stationCapacities))
        @assert(numAmbs <= sum(stationCapacities))
        remainingCapacity = copy(stationCapacities)
        unfilledStations = Set(findall(!iszero, remainingCapacity))
        priorityList = PriorityList(undef, numAmbs)
        priorityList[:] .= 0
        for i = 1:numAmbs
            j = rand(rng, unfilledStations) # station index
            priorityList[i] = j
            remainingCapacity[j] -= 1
            if remainingCapacity[j] == 0
                delete!(unfilledStations, j)
            end
        end
        return priorityList
    end
end
