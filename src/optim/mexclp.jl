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

# Maximum Expected Coverage Location Problem (MEXCLP)
# From paper: "A maximum expected covering location model: formulation, properties and heuristic solution"

"""
    solveMexclp!(sim::Simulation;
        numAmbs::Int = sim.numAmbs,
        busyFraction::Float = 0.5,
        demandWeights::Dict{Priority,Float} = Dict([p => 1.0 for p in priorities]),
        stationCapacities::Vector{Int} = [station.capacity for station in sim.stations],
        results::Dict = Dict())
Solves the Maximum Expected Coverage Location Problem (MEXCLP) for `sim` and returns the number of ambulances to assign to each station, also the converse is returned - a station index for each ambulance.
The problem assumes that all ambulances are equivalent.

# Keyword arguments
- `numAmbs` is the number of ambulances to solve for, must be >= 1
- `busyFraction` is the fraction of time that ambulances are busy; should be within [0,1] though this is not enforced
- `demandWeights` is the weight to apply to each demand priority for the objective function
- `stationCapacities` is the maximum number of ambulances that each station can hold
- `results` will store results of mexclp such as the objective value and decision variable values
"""
function solveMexclp!(sim::Simulation;
    numAmbs::Int=sim.numAmbs,
    busyFraction::Float=0.5,
    demandWeights::Dict{Priority,Float}=Dict([p => 1.0 for p in priorities]),
    stationCapacities::Vector{Int}=[station.capacity for station in sim.stations],
    results::Dict=Dict())

    @assert(numAmbs >= 1, "need at least 1 ambulance for mexclp")
    @assert(sim.travel.numSets == 1) # otherwise need to solve mexclp for each travel set?

    @assert(length(stationCapacities) == sim.numStations)
    @assert(all(stationCapacities .>= 0), "station capacities must be non-negative")
    stationCapacities = min.(stationCapacities, numAmbs) # reduce values where capacity > numAmbs
    @assert(sum(stationCapacities) >= numAmbs, "the total capacity for ambulances at stations is less than the number of ambulances")

    # initialise demand and demand coverage data if not already initialised
    sim.demand.initialised || initDemand!(sim)
    @assert(sim.demand.numSets == 1) # otherwise need to solve mexclp for each demand set?
    sim.demandCoverage.initialised || initDemandCoverage!(sim)

    # shorthand
    numStations = sim.numStations
    currentTime = sim.startTime

    # get demand point coverage data
    pointData = Dict{Vector{Int},Float}() # pointData[pointStations] = pointDemand
    for demandPriority in priorities
        if !haskey(demandWeights, demandPriority) || demandWeights[demandPriority] == 0
            continue
        end

        pointsCoverageMode = getPointsCoverageMode!(sim, demandPriority, currentTime)
        demandMode = getDemandMode!(sim.demand, demandPriority, currentTime)
        pointSetsDemands = getPointSetsDemands!(sim, demandPriority, currentTime; pointsCoverageMode=pointsCoverageMode) * demandMode.rasterMultiplier * demandWeights[demandPriority]

        for (i, stationSet) in enumerate(pointsCoverageMode.stationSets)
            get!(pointData, stationSet, 0.0)
            pointData[stationSet] += pointSetsDemands[i]
        end
    end
    pointStations = collect(keys(pointData))
    pointDemands = collect(values(pointData))
    # point j has demand pointDemands[j] and is covered by stations pointStations[j]
    numPoints = length(pointDemands) # shorthand

    # calculate benefit of covering each point with a kth ambulance
    marginalBenefit = (busyFraction .^ [0:numAmbs-1;]) * (1 - busyFraction) # point cover benefit values, for single demand
    pointCoverCountBenefit = [pointDemands[j] * marginalBenefit for j = 1:numPoints] # pointCoverCountBenefit[j][k] is the benefit of adding a kth ambulance to cover point j

    #######
    # IP

    # shorthand:
    a = numAmbs
    s = numStations
    p = numPoints

    model = Model()

    if !all(v -> issorted(v, rev=true), pointCoverCountBenefit)
        set_optimizer(model, optimizer_with_attributes(GLPK.Optimizer, "msg_lev" => GLPK.GLP_MSG_OFF)) # solve speed not tested

        @variables(model, begin
            (x[i=1:s] >= 0, Int) # x[i] = number of ambulances assigned to station i
            (y[j=1:p, k=1:a], Bin) # y[j,k] = true if point j is covered by at least k ambulances
        end)

        @constraints(model, begin
            (useAllAmbs, sum(x) == a)
            (pointCoverCount[j=1:p], sum(y[j, :]) == sum(x[pointStations[j]]))
            (pointCoverOrder[j=1:p, k=1:(a-1)], y[j, k] >= y[j, k+1])
            (stationMaxNumAmbs[i=1:s], x[i] <= stationCapacities[i])
        end)
    else
        # pointCoverCountBenefit[j] is non-increasing (because busyFraction >= 0), so:
        # - can have y variables as non-binary (but still constrained to be between 0 and 1)
        # - can leave out the constraint 'pointCoverOrder'
        # - could have x variables as non-integer, though if there are multiple solutions then x values might not be naturally integer, so will leave the integer constraint

        set_optimizer(model, optimizer_with_attributes(Cbc.Optimizer, "logLevel" => 0)) # solve speed not tested

        @variables(model, begin
            (x[i=1:s] >= 0, Int) # x[i] = number of ambulances assigned to station i
            (0 <= y[j=1:p, k=1:a] <= 1) # y[j,k] = true if point j is covered by at least k ambulances
        end)

        @constraints(model, begin
            (useAllAmbs, sum(x) == a)
            (pointCoverCount[j=1:p], sum(y[j, :]) == sum(x[pointStations[j]]))
            (stationMaxNumAmbs[i=1:s], x[i] <= stationCapacities[i])
        end)
    end

    @expressions(model, begin
        # (expectedPointCoverage[j=1:p], sum(y[j,k] * pointCoverCountBenefit[j][k] for k=1:a))
        (expectedCoverage, sum(y[j, k] * pointCoverCountBenefit[j][k] for j = 1:p, k = 1:a))
    end)

    # solve
    @objective(model, Max, expectedCoverage)
    @stdout_silent optimize!(model)
    @assert(termination_status(model) == MOI.OPTIMAL)
    xValue = JuMP.value.(x)
    yValue = JuMP.value.(y) # JuMP and LightXML both export value()
    objVal = JuMP.value(expectedCoverage)

    results[:x] = xValue
    results[:y] = yValue
    results[:objVal] = objVal

    # extract solution
    stationsNumAmbs = convert(Vector{Int}, round.(xValue)) # solution; stationsNumAmbs[i] gives number of ambulances to allocate to station i
    pointSlotCover = convert(Array{Bool,2}, round.(yValue)) # pointSlotCover[j,k] = true if point j is covered by at least k ambulances

    if checkMode
        # check constraints
        @assert(sum(stationsNumAmbs) == a) # useAllAmbs constraint
        @assert(all(j -> sum(pointSlotCover[j, :]) == sum(stationsNumAmbs[pointStations[j]]), 1:p)) # pointCoverCount constraint
        @assert(all(i -> 0 <= stationsNumAmbs[i] <= stationCapacities[i], 1:s)) # stationMaxNumAmbs constraint
        @assert(all(j -> issorted(pointSlotCover[j], rev=true), 1:p)) # pointCoverOrder constraint
    end

    # convert to a deployment
    deployment = stationsNumAmbsToDeployment(stationsNumAmbs)

    return stationsNumAmbs, deployment # also return model?
end
