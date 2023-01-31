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

# Maximal Covering Location Problem (MCLP)
# From paper: "The maximal covering location problem"

"""
    solveMclp!(sim::Simulation;
        numAmbs::Int = sim.numAmbs,
        demandWeights::Dict{Priority,Float} = Dict([p => 1.0 for p in priorities]),
        results::Dict = Dict())
Solves the Maximal Covering Location Problem (MCLP) to locate `numAmbs` ambulances at stations for `sim` and returns the number of ambulances to assign to each station, also the converse is returned - a station index for each ambulance.
The problem assumes that all ambulances are equivalent.

# Keyword arguments
- `numAmbs` is the number of ambulances to solve for, should be in `1:sim.numStations`
- `demandWeights` is the weight to apply to each demand priority for the objective function
- `results` will store results of mclp such as the objective value and decision variable values
"""
function solveMclp!(sim::Simulation;
    numAmbs::Int=sim.numAmbs,
    demandWeights::Dict{Priority,Float}=Dict([p => 1.0 for p in priorities]),
    results::Dict=Dict())

    @assert(1 <= numAmbs <= sim.numStations, "mclp requires that numAmbs be in 1:sim.numStations")
    @assert(sim.travel.numSets == 1) # otherwise need to solve mclp for each travel set?

    # initialise demand and demand coverage data if not already initialised
    sim.demand.initialised || initDemand!(sim)
    @assert(sim.demand.numSets == 1) # otherwise need to solve mclp for each demand set?
    sim.demandCoverage.initialised || initDemandCoverage!(sim)

    # get demand point coverage data
    currentTime = sim.startTime # shorthand
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

    stationsNumAmbs = solveMclp(numAmbs, pointDemands, pointStations; results=results)
    deployment = stationsNumAmbsToDeployment(stationsNumAmbs) # convert to a deployment

    return stationsNumAmbs, deployment
end

"""
    solveMclp(n::Int, pointDemands::Vector{Float}, coverMatrix::Array{Bool,2}; results::Dict = Dict())
Solves the Maximal Covering Location Problem (MCLP) to locate `n` facilities, where `pointDemands[j]` is the demand at point `j` and `coverMatrix[i,j] = true` if facility location `i` can cover point `j`.

`results` will store results of mclp such as the objective value and decision variable values (though the decision variables may have changed from the given problem as duplicate columns in `coverMatrix` are removed).
"""
function solveMclp(n::Int, pointDemands::Vector{Float}, coverMatrix::Array{Bool,2}; results::Dict=Dict())
    l = size(coverMatrix, 1) # number of potential facility locations
    @assert(all(x -> x >= 0, pointDemands))
    @assert(0 <= n <= l, "mclp requires that n be between 0 and the number of potential facility locations ($l).")
    @assert(length(pointDemands) == size(coverMatrix, 2))

    # combine points that are covered by the same facilities
    pointData = Dict{Vector{Bool},Float}() # pointData[coverMatrixCol] = pointSetDemand
    for (j, pointDemand) in enumerate(pointDemands)
        coverMatrixCol = coverMatrix[:, j]
        get!(pointData, coverMatrixCol, 0.0)
        pointData[coverMatrixCol] += pointDemand
    end

    # Could also combine any facility locations that cover the same points,
    # or remove facility locations with coverage dominated by another location,
    # but I'll leave that as it is unlikely for my use cases.

    # # remove points that can be covered by more facility locations than the number of facilities that will not be opened,
    # # as it is guaranteed that these points will be covered
    # objValDelta = 0.0
    # for (coverMatrixCol, pointDemand) in pointData
    # if sum(coverMatrixCol) > l - n
    # delete!(pointData, coverMatrixCol)
    # objValDelta += pointDemand
    # end
    # end

    pointSetDemands = collect(values(pointData))
    pointSetFacilities = findall.(collect(keys(pointData))) # pointSetFacilities[i] gives indices of facility locations that cover pointSetDemands[i]

    # # remove point set that cannot be covered
    # i = findfirst(isempty, pointSetFacilities)
    # if i != nothing
    # deleteat!(pointSetFacilities, i)
    # deleteat!(pointSetDemands, i)
    # end

    solution = solveMclp(n, pointSetDemands, pointSetFacilities, results=results)
    solution = vcat(solution, zeros(Int, size(coverMatrix, 1) - length(solution))) # right-pad with zeros if needed

    # # change results dict to be for values given from function args instead of modified problem (that will be same size as given, or smaller)
    # results[:objVal] += objValDelta
    # # results[:x] should be unchanged, unless bottom row(s) of coverMatrix are all false
    # # results[:y] will not be for given pointDemands unless all columns of coverMatrix are unique

    return solution
end

"""
    solveMclp(n::Int, pointDemands::Vector{Float}, pointFacilities::Vector{Vector{Int}}; results::Dict = Dict())
Solves the Maximal Covering Location Problem (MCLP) to locate `n` facilities, where `pointDemands[i]` is the demand at point `i` and `pointFacilities[i]` are the indices of facility locations that cover point `i`.

`results` will store results of mclp such as the objective value and decision variable values.
"""
function solveMclp(n::Int, pointDemands::Vector{Float}, pointFacilities::Vector{Vector{Int}}; results::Dict=Dict())
    p = length(pointDemands) # number of demand points
    @assert(all(x -> x >= 0, pointDemands))
    l = maximum([isempty(v) ? 0 : maximum(v) for v in pointFacilities]) # number of potential facility locations
    @assert(0 <= n <= l, "mclp requires that n be between 0 and the number of potential facility locations ($l).")
    @assert(length(pointFacilities) == p)

    model = Model()

    # use cbc solver; solve speed not tested
    set_optimizer(model, with_optimizer(Cbc.Optimizer, logLevel=0))

    @variables(model, begin
        (x[i=1:l], Bin) # x[i] = true if facility location i is used, false otherwise
        (0 <= y[j=1:p] <= 1) # y[j] = true if point j is covered
    end)

    @constraints(model, begin
        (useAllFacilities, sum(x) == n)
        (pointCoverage[j=1:p], y[j] <= sum(x[pointFacilities[j]]))
    end)

    @expressions(model, begin
        (coverage, sum(y .* pointDemands))
    end)

    # solve
    @objective(model, Max, coverage)
    optimize!(model)
    @assert(termination_status(model) == MOI.OPTIMAL)
    xValue = JuMP.value.(x)
    yValue = JuMP.value.(y) # JuMP and LightXML both export value()
    objVal = JuMP.value(coverage)

    results[:x] = xValue
    results[:y] = yValue
    results[:objVal] = objVal

    # extract solution
    solution = convert(Vector{Int}, round.(xValue))

    if checkMode
        # check constraints
        @assert(sum(solution) == n) # useAllFacilities constraint
    end

    return solution
end
