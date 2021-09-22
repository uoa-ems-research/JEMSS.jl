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

# multi-compliance table for move up

# initialise data relevant to move up
function initMultiCompTable!(sim::Simulation, multiCompTable::MultiCompTable)
	# shorthand names:
	@unpack numAmbs, numStations = sim
	mctd = sim.moveUpData.multiCompTableData
	
	# set multi-compliance table
	checkMultiCompTable(multiCompTable, sim) # not checking station capacities
	mctd.multiCompTable = multiCompTable
end
function initMultiCompTable!(sim::Simulation, multiCompTableFilename::String)
	multiCompTable = readMultiCompTableFile(multiCompTableFilename)
	initMultiCompTable!(sim, multiCompTable)
end

function multiCompTableMoveUp(sim::Simulation)
	@assert(sim.moveUpData.useMoveUp)
	
	# shorthand
	@unpack ambulances, stations, numAmbs, numStations = sim
	mctd = sim.moveUpData.multiCompTableData
	multiCompTable = mctd.multiCompTable
	
	# get movable ambulances (movableAmbs)
	ambMovable = isAmbMovable.(ambulances)
	movableAmbs = ambulances[ambMovable]
	numMovableAmbs = length(movableAmbs)
	
	if numMovableAmbs == 0
		return moveUpNull()
	end
	
	# calculate travel time for each movable ambulance
	ambToStationTimes = zeros(Float, numMovableAmbs, numStations) # ambToStationTimes[i,j] = time for movable ambulance i to travel to station j
	for i = 1:numMovableAmbs
		ambToStationTimes[i,:] = ambMoveUpTravelTimes!(sim, movableAmbs[i])
	end
	
	# for all applicable comp table states, solve assignment problem and calculate total travel time
	ambStationIndicesList = []
	totalTravelTimes = []
	compTableStates = multiCompTable[numMovableAmbs]
	for compTableState in compTableStates
		stationSlots = vcat([[stationIndex for j = 1:stationSlotCount] for (stationIndex, stationSlotCount) in enumerate(compTableState)]...)
		results = Dict()
		ambStationIndices = solveCompTableAssignmentProblem(stationSlots, ambToStationTimes, results = results)
		push!(ambStationIndicesList, ambStationIndices)
		push!(totalTravelTimes, results[:objVal])
	end
	
	# use comp table state with smallest total travel time
	i = argmin(totalTravelTimes)
	ambStationIndices = ambStationIndicesList[i]
	compTableState = compTableStates[i]
	
	ambStations = stations[ambStationIndices]
	
	if checkMode
		# check that comp table state is followed
		compTableStateCopy = deepcopy(compTableState)
		for j in ambStationIndices
			compTableStateCopy[j] -= 1
		end
		@assert(iszero(compTableStateCopy))
	end
	
	return movableAmbs, ambStations
end

# check that multi-compliance table is valid
function checkMultiCompTable(multiCompTable::MultiCompTable;
	numAmbs::Int = nullIndex, numStations::Int = nullIndex,
	stationCapacities::Union{Vector{Int},Nothing} = nothing)
	
	@assert(numAmbs == nullIndex || length(multiCompTable) == numAmbs)
	if numStations == nullIndex
		numStations = multiCompTable[1][1]
	end
	for i = 1:numAmbs
		states = multiCompTable[i]
		@assert(length(states) >= 1)
		for state in states
			# state[k] gives number of ambulances to assign to station k
			@assert(length(state) == numStations)
			@assert(all(x -> x >= 0, state))
			@assert(sum(state) == i) # match number of ambs
			if stationCapacities !== nothing
				@assert(all(state .<= stationCapacities))
			end
		end
	end
	return true
end
function checkMultiCompTable(multiCompTable::MultiCompTable, sim::Simulation)
	return checkMultiCompTable(multiCompTable; numAmbs = sim.numAmbs, numStations = sim.numStations, stationCapacities = [s.capacity for s in sim.stations])
end
