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
function initCompTable!(sim::Simulation, compTableFilename::String)
	# shorthand names:
	numAmbs = sim.numAmbs
	numStations = sim.numStations
	ctd = sim.moveUpData.compTableData
	
	# read compliance table
	ctd.compTable = readCompTableFile(compTableFilename)
	@assert(size(ctd.compTable) == (numAmbs, numStations))
	
	# compliance table should not violate station capacities
	for j = 1:numStations
		@assert(all(ctd.compTable[:,j] .<= sim.stations[j].capacity))
	end
	
	ctd.compTableStationSlots = [vcat([[stationIndex for j = 1:stationSlotCount] for (stationIndex, stationSlotCount) in enumerate(ctd.compTable[i,:])]...) for i = 1:numAmbs] # compTableStationSlots[i] gives the indices of stations as many times as the number of ambulances required at the station for compTable[i,:]
	
	# set array sizes
	ctd.ambMovable = Vector{Bool}(numAmbs)
end

function compTableMoveUp(sim::Simulation)
	@assert(sim.moveUpData.useMoveUp)
	
	# shorthand:
	ambulances = sim.ambulances
	stations = sim.stations
	numAmbs = sim.numAmbs
	numStations = sim.numStations
	ctd = sim.moveUpData.compTableData
	compTable = ctd.compTable
	ambMovable = ctd.ambMovable
	
	# get movable ambulances (movableAmbs)
	for i = 1:numAmbs
		ambMovable[i] = isAmbAvailableForMoveUp(ambulances[i])
	end
	movableAmbs = ambulances[ambMovable]
	numMovableAmbs = length(movableAmbs)
	
	if numMovableAmbs == 0
		return [], []
	end
	
	# calculate travel time for each available ambulance to reach every station that requires >= 1 amb
	ambToStationTimes = zeros(Float, numMovableAmbs, numStations) # ambToStationTimes[i,j] = time for movable ambulance i to travel to station j
	moveUpStations = find(compTable[numMovableAmbs,:]) # indices of stations that require >= 1 amb
	for i = 1:numMovableAmbs
		ambToStationTimes[i,moveUpStations] = ambMoveUpTravelTimes!(sim, movableAmbs[i]; stations = stations[moveUpStations])
	end
	
	# solve as an integer program
	# ambStationIndices = solveCompTableIP(compTable[numMovableAmbs,:], ambToStationTimes)
	
	# solve as an assignment problem
	ambStationIndices = solveCompTableAssignmentProblem(ctd.compTableStationSlots[numMovableAmbs], ambToStationTimes)
	
	ambStations = stations[ambStationIndices]
	
	if checkMode
		# check that compTable is followed
		compTableRow = compTable[numMovableAmbs,:]
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
	
	model = Model()
	
	setsolver(model, GLPKSolverMIP(presolve=true))
	
	@variables(model, begin
		(x[i=1:a,j=1:s], Bin) # x[i,j] = 1 if ambulance: movableAmbs[i] should be moved to station: stations[j]
		# maxTravelCost # if minimising max travel cost of all ambs, instead of total
	end)
	
	@constraints(model, begin
		(followCompTable[j=1:s], sum(x[i,j] for i=1:a) == compTableRow[j])
		(ambAtOneLocation[i=1:a], sum(x[i,j] for j=1:s) == 1) # each ambulance must be assigned to one station
		# maxTravelCostBound[i=1:a], sum(x[i,j] * ambToStationCost[i,j] for j=1:s) <= maxTravelCost # max travel cost of all ambulances
	end)
	
	@expressions(model, begin
		totalAmbTravelCost, sum(x[i,j] * ambToStationCost[i,j] for i=1:a, j=1:s)
	end)
	
	@objective(model, :Min, totalAmbTravelCost)
	# @objective(model, :Min, maxTravelCost)
	
	# # testing: giving back fake results, for testing runtime without solving IP model
	# if true
		# return moveUpNull()
	# end
	
	solve(model)
	
	# extract solution
	sol = convert(Array{Bool,2}, round.(getvalue(x)))
	(ambIndices, stationIndices) = findn(sol)
	ambStationIndices = zeros(Int,numMovableAmbs) # ambStationIndices[i] gives index of the station that movable ambulance i should be assigned to
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
function solveCompTableAssignmentProblem(stationSlots::Vector{Int}, ambToStationCost::Array{Float,2})
	# stationSlots gives the indices of stations as many times as the number of ambulances required at the station
	# ambToStationCost[i,j] is the cost of assigning movable ambulance i to station j
	# formulated as an assignment problem and solved with the Hungarian algorithm
	
	(numMovableAmbs, numStations) = size(ambToStationCost)
	@assert(length(stationSlots) == numMovableAmbs)
	@assert(all(stationIndex -> 1 <= stationIndex <= numStations, stationSlots))
	
	# formulate assignment problem and solve
	# need to create as many copies of station j according to values in stationSlots
	weights = ambToStationCost[:,stationSlots] # size(weights) is (numMovableAmbs, numMovableAmbs)
	matching = Hungarian.munkres(weights) # returns a sparse matrix
	# (assignment, cost) = hungarian(weights) # slightly slower than Hungarian.munkres(weights); it allows for dummy nodes but I don't need this functionality
	
	# extract solution
	(ambIndices, stationSlotIndices) = findn(matching .== Hungarian.STAR)
	ambStationIndices = zeros(Int,numMovableAmbs) # ambStationIndices[i] gives index of the station that movable ambulance i should be assigned to
	ambStationIndices[ambIndices] = stationSlots[stationSlotIndices]
	
	# if checkMode
		# check that compliance table is followed
		# @assert(Set(ambStationIndices) == Set(stationSlots)) # this may be slow
	# end
	
	return ambStationIndices
end
