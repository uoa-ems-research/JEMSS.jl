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

# simulation replications

# mutates: sim, callSets
function setSimReps!(sim::Simulation, callSets::Vector{Vector{Call}})
	sim.reps = []
	for calls in callSets
		initCalls!(sim, calls)
		rep = Simulation()
		rep.calls = calls
		rep.numCalls = length(calls)
		push!(sim.reps, rep)
	end
	sim.numReps = length(sim.reps)
	if isdefined(sim, :backup) sim.backup.numReps = sim.numReps end
	
	sim.writeOutput = false # have not yet handled writing output files for each replication
end

# run single simulation replication
function simulateRep!(sim::Simulation, rep::Simulation)
	sim.calls = Call[] # so calls will not be reset
	reset!(sim)
	setSimCalls!(sim, rep.calls)
	simulate!(sim)
	for fname in (:startTime, :time, :endTime, :stats, :used, :complete)
		setfield!(rep, fname, getfield(sim, fname))
	end
end
simulateRep!(sim::Simulation, repIndex::Int) = simulateRep!(sim, sim.reps[repIndex])

# run simulation replications
function simulateReps!(sim::Simulation; repIndices::Vector{Int} = [1:sim.numReps;], doPrint::Bool = false)
	@assert(all(in(1:sim.numReps), repIndices))
	reps = sim.reps[repIndices] # shorthand
	numReps = length(reps)
	for (i,rep) in enumerate(reps)
		doPrint && print("\rSimulating replication $i of $numReps.")
		simulateRep!(sim, rep)
	end
	doPrint && println()
end
