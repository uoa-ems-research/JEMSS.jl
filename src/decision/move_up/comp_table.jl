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

# compliance table for move up

# initialise data relevant to move up
function initCompTable!(sim::Simulation, compTable::CompTable)
    # shorthand names:
    @unpack numAmbs, numStations = sim
    ctd = sim.moveUpData.compTableData

    # set compliance table
    checkCompTable(compTable, sim) # not checking station capacities
    ctd.compTable = compTable

    ctd.compTableStationSlots = [vcat([[stationIndex for j = 1:stationSlotCount] for (stationIndex, stationSlotCount) in enumerate(ctd.compTable[i, :])]...) for i = 1:numAmbs] # compTableStationSlots[i] gives the indices of stations as many times as the number of ambulances required at the station for compTable[i,:]

    # set array sizes
    ctd.ambMovable = Vector{Bool}(undef, numAmbs)
end
function initCompTable!(sim::Simulation, nestedCompTable::NestedCompTable)
    compTable = unnestCompTable(nestedCompTable, sim.numStations)
    initCompTable!(sim, compTable)
end
function initCompTable!(sim::Simulation, compTableFilename::String)
    compTable = readCompTableFile(compTableFilename)
    initCompTable!(sim, compTable)
end

function compTableMoveUp(sim::Simulation)
    @assert(sim.moveUpData.useMoveUp)

    # shorthand:
    @unpack ambulances, stations, numAmbs, numStations = sim
    ctd = sim.moveUpData.compTableData
    @unpack compTable, ambMovable = ctd

    # get movable ambulances (movableAmbs)
    for i = 1:numAmbs
        ambMovable[i] = isAmbMovable(ambulances[i])
    end
    movableAmbs = ambulances[ambMovable]
    numMovableAmbs = length(movableAmbs)

    if numMovableAmbs == 0
        return moveUpNull()
    end

    # calculate travel time for each movable ambulance to reach every station that requires >= 1 amb
    ambToStationTimes = zeros(Float, numMovableAmbs, numStations) # ambToStationTimes[i,j] = time for movable ambulance i to travel to station j
    moveUpStations = findall(!iszero, compTable[numMovableAmbs, :]) # indices of stations that require >= 1 amb
    for i = 1:numMovableAmbs
        ambToStationTimes[i, moveUpStations] = ambMoveUpTravelTimes!(sim, movableAmbs[i]; stations=stations[moveUpStations])
    end

    # solve as an integer program
    # ambStationIndices = solveCompTableIP(compTable[numMovableAmbs,:], ambToStationTimes)

    # solve as an assignment problem
    ambStationIndices = solveCompTableAssignmentProblem(ctd.compTableStationSlots[numMovableAmbs], ambToStationTimes)

    ambStations = stations[ambStationIndices]

    if checkMode
        # check that compTable is followed
        compTableRow = compTable[numMovableAmbs, :]
        for j in ambStationIndices
            compTableRow[j] -= 1
        end
        @assert(iszero(compTableRow))
    end

    return movableAmbs, ambStations
end

function solveCompTableIP(compTableRow::Vector{Int}, ambToStationCost::Array{Float,2})
    # solve compliance table problem as an IP, minimising total cost
    # compTableRow[j] gives the number of ambulances to allocate to station j
    # ambToStationCost[i,j] is the cost of assigning movable ambulance i to station j

    # shorthand:
    a = numMovableAmbs = sum(compTableRow)
    s = numStations = length(compTableRow)

    # using JuMP
    model = Model()
    set_optimizer(model, with_optimizer(GLPK.Optimizer)) # solve speed not tested
    # set_optimizer(model, with_optimizer(Cbc.Optimizer)) # not sure how to mute output from cbc

    @variables(model, begin
        (x[i=1:a, j=1:s], Bin) # x[i,j] = 1 if ambulance: movableAmbs[i] should be moved to station: stations[j]
        # maxTravelCost # if minimising max travel cost of all ambs, instead of total
    end)

    @constraints(model, begin
        (followCompTable[j=1:s], sum(x[i, j] for i = 1:a) == compTableRow[j])
        (ambAtOneLocation[i=1:a], sum(x[i, j] for j = 1:s) == 1) # each ambulance must be assigned to one station
        # maxTravelCostBound[i=1:a], sum(x[i,j] * ambToStationCost[i,j] for j=1:s) <= maxTravelCost # max travel cost of all ambulances
    end)

    @expressions(model, begin
        totalAmbTravelCost, sum(x[i, j] * ambToStationCost[i, j] for i = 1:a, j = 1:s)
    end)

    # # testing: giving back fake results, for testing runtime without solving IP model
    # if true
    # return moveUpNull()
    # end

    # solve
    @objective(model, Min, totalAmbTravelCost) # or maxTravelCost
    optimize!(model)
    @assert(termination_status(model) == MOI.OPTIMAL)
    xValue = JuMP.value.(x) # JuMP and LightXML both export value()

    # extract solution
    sol = convert(Array{Bool,2}, round.(xValue))
    I = findall(sol)
    (ambIndices, stationIndices) = (getindex.(I, 1), getindex.(I, 2))
    ambStationIndices = zeros(Int, numMovableAmbs) # ambStationIndices[i] gives index of the station that movable ambulance i should be assigned to
    ambStationIndices[ambIndices] = stationIndices

    # if checkMode
    # # check that followCompTable constraint was met
    # for j = 1:numStations
    # @assert(sum(sol[:,j]) == compTableRow[j])
    # end
    # # check that ambAtOneLocation constraint was met
    # for i = 1:numMovableAmbs
    # @assert(sum(sol[i,:]) == 1)
    # end
    # end

    return ambStationIndices
end

# solve the compliance table problem with the Hungarian algorithm (for the assignment problem),
# this can only be used if the objective is to minimise the sum of individual ambulance redeployment costs
function solveCompTableAssignmentProblem(stationSlots::Vector{Int}, ambToStationCost::Array{Float,2}; results::Dict=Dict())
    # stationSlots gives the indices of stations as many times as the number of ambulances required at the station
    # ambToStationCost[i,j] is the cost of assigning movable ambulance i to station j
    # formulated as an assignment problem and solved with the Hungarian algorithm

    (numMovableAmbs, numStations) = size(ambToStationCost)
    @assert(length(stationSlots) == numMovableAmbs)
    @assert(all(stationIndex -> 1 <= stationIndex <= numStations, stationSlots))

    # formulate assignment problem and solve
    # need to create as many copies of station j according to values in stationSlots
    weights = ambToStationCost[:, stationSlots] # size(weights) is (numMovableAmbs, numMovableAmbs)
    matching = Hungarian.munkres(weights) # returns a sparse matrix
    # (assignment, cost) = hungarian(weights) # slightly slower than Hungarian.munkres(weights); it allows for dummy nodes but I don't need this functionality

    # extract solution
    I = findall(matching .== Hungarian.STAR)
    (ambIndices, stationSlotIndices) = (getindex.(I, 1), getindex.(I, 2))
    ambStationIndices = zeros(Int, numMovableAmbs) # ambStationIndices[i] gives index of the station that movable ambulance i should be assigned to
    ambStationIndices[ambIndices] = stationSlots[stationSlotIndices]

    # if checkMode
    # check that compliance table is followed
    # @assert(Set(ambStationIndices) == Set(stationSlots)) # this may be slow
    # end

    results[:objVal] = sum(weights[I])

    return ambStationIndices
end

# check that compliance table is valid
function checkCompTable(compTable::CompTable;
    numAmbs::Int=nullIndex, numStations::Int=nullIndex,
    stationCapacities::Union{Vector{Int},Nothing}=nothing)

    (m, n) = size(compTable)
    @assert(all(v -> v >= 0, compTable)) # need non-negative integers
    for i = 1:m
        @assert(sum(compTable[i, :]) == i) # row sum should match row index (= number of movable ambulances)
    end
    @assert(numAmbs == nullIndex || numAmbs == m)
    @assert(numStations == nullIndex || numStations == n)
    if stationCapacities !== nothing
        @assert(length(stationCapacities) == n)
        for j = 1:n
            @assert(all(v -> v <= stationCapacities[j], compTable[:, j]))
            # @assert(all(compTable[:,j] .<= stationCapacities[j]))
        end
    end
    return true
end
function checkCompTable(compTable::CompTable, sim::Simulation)
    return checkCompTable(compTable; numAmbs=sim.numAmbs, numStations=sim.numStations, stationCapacities=[s.capacity for s in sim.stations])
end

function isCompTableNested(compTable::CompTable)
    checkCompTable(compTable)
    for i = 1:size(compTable)-1
        if count(!iszero, compTable[i+1, :] - compTable[i, :]) != 1 # should only have one change between rows
            return false
        end
    end
    return true
end

# take a compliance table and return the nested form, if possible
function nestCompTable(compTable::CompTable)::NestedCompTable
    @assert(isCompTableNested(compTable)) # check if nesting is possible
    (m, n) = size(compTable)
    nestedCompTable = NestedCompTable(undef, m)
    nestedCompTable[1] = findfirst(!iszero, compTable[1, :])
    for i = 2:m
        nestedCompTable[i] = findfirst(!iszero, compTable[i, :] - compTable[i-1, :])
    end
    return nestedCompTable
end

# convert nested comp table to normal compliance table (list to array)
function unnestCompTable(nestedCompTable::NestedCompTable, numStations::Int)::CompTable
    m = length(nestedCompTable)
    n = numStations
    compTable = CompTable(undef, m, n)
    compTable[:] .= 0
    for i = 1:m
        compTable[i:m, nestedCompTable[i]] .+= 1
    end
    return compTable
end

# returns a randomly generated nested comp table
function makeRandNestedCompTable(numAmbs::Int, numStations::Int;
    stationCapacities::Union{Vector{Int},Nothing}=nothing,
    rng::AbstractRNG=GLOBAL_RNG)::NestedCompTable

    @assert(numStations > 0)
    if stationCapacities === nothing
        return rand(rng, 1:numStations, numAmbs)
    else
        @assert(numStations == length(stationCapacities))
        @assert(numAmbs <= sum(stationCapacities))
        remainingCapacity = copy(stationCapacities)
        unfilledStations = Set(findall(!iszero, remainingCapacity))
        nestedCompTable = NestedCompTable(undef, numAmbs)
        nestedCompTable[:] .= 0
        for i = 1:numAmbs
            j = rand(rng, unfilledStations) # station index
            nestedCompTable[i] = j
            remainingCapacity[j] -= 1
            if remainingCapacity[j] == 0
                delete!(unfilledStations, j)
            end
        end
        return nestedCompTable
    end
end
