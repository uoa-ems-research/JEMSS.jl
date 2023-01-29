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

# Maximal Coverage Location Problem applied to creating a nested compliance table

"""
	solveNestedMclp!(sim::Simulation, numFreeAmbsWeights::Vector{Float};
		demandWeights::Dict{Priority,Float} = Dict([p => 1.0 for p in priorities]),
		results::Dict{Any,Any} = Dict())
Solves combined Maximal Coverage Location Problems (MCLP) to create a nested compliance table.
`numFreeAmbsWeights[k]` gives the weight to apply to the objective function for MCLP with `k` ambulances; this can be used to reflect the probability of having `k` ambulances free at any moment.
Ambulance travel time between states (i.e. rows) of the compliance table is not accounted for.

# Keyword arguments
- `demandWeights` is the weight to apply to each demand priority for the objective function
- `results` will store results such as the objective value and decision variable values
"""
function solveNestedMclp!(sim::Simulation, numFreeAmbsWeights::Vector{Float};
    demandWeights::Dict{Priority,Float}=Dict([p => 1.0 for p in priorities]),
    results::Dict{Any,Any}=Dict())

    @assert(all(numFreeAmbsWeights .>= 0))

    @assert(sim.demand.numSets == 1) # otherwise need to solve mclp for each demand set?
    @assert(sim.travel.numSets == 1) # otherwise need to solve mclp for each travel set?

    # shorthand
    maxNumAmbs = sim.numStations # any more ambs than this and MCLP solutions are equivalent
    currentTime = sim.startTime

    if length(numFreeAmbsWeights) < maxNumAmbs
        numFreeAmbsWeights = vcat(numFreeAmbsWeights, zeros(Float, maxNumAmbs - length(numFreeAmbsWeights)))
    else
        # change numFreeAmbsWeights to reflect that MCLP solution does not change when the number of free ambs >= num stations
        numFreeAmbsWeights = vcat(numFreeAmbsWeights[1:maxNumAmbs-1], sum(numFreeAmbsWeights[maxNumAmbs:end]))
    end
    @assert(length(numFreeAmbsWeights) == maxNumAmbs)

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

    nestedCompTable = solveNestedMclp(pointDemands, pointStations, numFreeAmbsWeights, results=results)

    return nestedCompTable
end

"""
	solveNestedMclp(pointDemands::Vector{Float}, coverMatrix::Array{Bool,2}, weights::Vector{Float}; results::Dict{Any,Any} = Dict())
Solves combined Maximal Coverage Location Problems (MCLP) to create a nested compliance table.
`pointDemands[i]` is the demand at point `i` and `coverMatrix[i,j] = true` if facility location `i` can cover point `j`.
`weights[k]` gives the weight to apply to the objective function for MCLP with `k` facilities.
`results` will store results such as the objective value and decision variable values.
"""
function solveNestedMclp(pointDemands::Vector{Float}, coverMatrix::Array{Bool,2}, weights::Vector{Float}; results::Dict{Any,Any}=Dict())
    l = size(coverMatrix, 1) # number of potential facility locations
    @assert(all(x -> x >= 0, pointDemands))
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

    pointSetDemands = collect(values(pointData))
    pointSetFacilities = findall.(collect(keys(pointData))) # pointSetFacilities[i] gives indices of facility locations that cover pointSetDemands[i]

    solution = solveNestedMclp(pointSetDemands, pointSetFacilities, weights, results=results)
    solution = vcat(solution, zeros(Int, size(coverMatrix, 1) - length(solution))) # right-pad with zeros if needed

    # change results dict to be for values given from function args instead of modified problem (that will be same size as given, or smaller)
    # results[:objVal] should be unchanged
    # results[:x] should be unchanged, unless bottom row(s) of coverMatrix are all false
    # results[:y] will not be for given pointDemands unless all columns of coverMatrix are unique

    return solution
end

"""
	solveNestedMclp(pointDemands::Vector{Float}, pointFacilities::Vector{Vector{Int}}, weights::Vector{Float}; results::Dict{Any,Any} = Dict())
Solves combined Maximal Coverage Location Problems (MCLP) to create a nested compliance table.
`pointDemands[i]` is the demand at point `i` and `pointFacilities[i]` are the indices of facility locations that cover point `i`.
`weights[k]` gives the weight to apply to the objective function for MCLP with `k` facilities.
`results` will store results such as the objective value and decision variable values.
"""
function solveNestedMclp(pointDemands::Vector{Float}, pointFacilities::Vector{Vector{Int}}, weights::Vector{Float}; results::Dict{Any,Any}=Dict())
    @assert(all(x -> x >= 0, pointDemands))
    p = length(pointDemands) # number of demand points
    @assert(length(pointFacilities) == p)
    @assert(all(weights .>= 0))
    f = length(weights) # number of potential facility locations
    @assert(all(x -> all(y -> in(y, 1:f), x), pointFacilities))

    #######
    # IP

    model = Model()

    # use cbc solver; solver speed not tested
    set_optimizer(model, with_optimizer(Cbc.Optimizer, logLevel=0))

    @variables(model, begin
        (x[i=1:f, k=1:f], Bin) # x[i,k] = 1 if facility location i is used when there are k facilities, 0 otherwise
        (0 <= y[j=1:p, k=1:f] <= 1) # y[j,k] = 1 if point j is covered when there are k facilities, 0 otherwise
    end)

    @constraints(model, begin
        (facilityCounts[k=1:f], sum(x[:, k]) == k) # can only select k facilities
        (pointCoverage[j=1:p, k=1:f], y[j, k] <= sum(x[pointFacilities[j], k])) # demand point is covered if at least one facility location that can cover it is used
        (nesting[i=1:f, k=1:f-1], x[i, k] <= x[i, k+1]) # if facility location i is used when there are k facilities, then it must also be used when there are k+1 facilities
    end)

    @expressions(model, begin
        (numFacilitiesCoverage[k=1:f], sum(y[j, k] * pointDemands[j] for j = 1:p))
        (weightedCoverage, sum(numFacilitiesCoverage[k] * weights[k] for k = 1:f))
    end)

    # solve
    @objective(model, Max, weightedCoverage)
    optimize!(model)
    @assert(termination_status(model) == MOI.OPTIMAL)
    xValue = JuMP.value.(x)
    yValue = JuMP.value.(y) # JuMP and LightXML both export value()
    objVal = JuMP.value(weightedCoverage)
    results[:numFacilitiesCoverage] = JuMP.value.(numFacilitiesCoverage)

    results[:objVal] = objVal
    results[:x] = xValue
    results[:y] = yValue

    # extract solution
    compTable = convert(Array{Int,2}, round.(xValue)') # solution; compTable[k,i] is 1 if facility location i is used when there are k facilities, 0 otherwise

    if checkMode
        # check constraints
        for i = 1:f
            @assert(sum(compTable[i, :]) == i) # row i of compTable should use i facilities
        end
    end

    return nestCompTable(compTable)
    # return compTable # alternatively, could return nestedCompTable, but it is easier to convert from compTable to nestedCompTable than convert the other way
end
