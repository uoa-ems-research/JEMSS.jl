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

# set simulation replications, one rep per call set
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
    if isdefined(sim, :backup)
        sim.backup.numReps = sim.numReps
    end

    # have not yet handled writing output files for each replication
    sim.writeOutput = false
    if isdefined(sim, :backup)
        sim.backup.writeOutput = false
    end
end

# run single simulation replication on `sim`
# need this function in case rep.isRunnable == false
# mutates: sim, rep
function simulateRep!(sim::Simulation, rep::Simulation)
    sim.calls = Call[] # so calls will not be reset
    reset!(sim)
    setSimCalls!(sim, rep.calls)
    simulate!(sim)
    for fname in (:startTime, :time, :endTime, :stats, :used, :complete)
        setfield!(rep, fname, getfield(sim, fname))
    end
end

# run simulation replications on `sim`
# need this function in case rep.isRunnable == false
# mutates: sim, reps
function simulateReps!(sim::Simulation, reps::Vector{Simulation}=sim.reps; doPrint::Bool=false)
    numReps = length(reps)
    for (i, rep) in enumerate(reps)
        doPrint && print("\rSimulating replication $i of $numReps.")
        simulateRep!(sim, rep)
    end
    doPrint && println()
end

# make sim replications runnable, e.g. in order to run simulate!(sim.reps[1])
# for experiments, any changes made to `sim` may require the same change to be made to each replication
# mutates: sim.reps
function makeRepsRunnable!(sim::Simulation)
    @assert(sim.initialised && !sim.used)
    for rep in sim.reps
        fnamesDontCopy = (:backup, :net, :travel, :grid, :resim, :calls, :demand, :demandCoverage, :reps)
        for fname in setdiff(fieldnames(Simulation), fnamesDontCopy)
            setfield!(rep, fname, deepcopy(getfield(sim, fname)))
        end

        # have some fields in rep point to same data as sim
        for fname in (:net, :grid, :demandCoverage)
            setfield!(rep, fname, getfield(sim, fname))
        end

        # travel and demand types have some fields that change value while simulating
        # in particular: sim.travel.recentSetsStartTimesIndex (Int) and sim.demand.recentSetsStartTimesIndex (Int)
        for fname in (:travel, :demand)
            dst = getfield(rep, fname)
            src = getfield(sim, fname)
            for fname2 in fieldnames(typeof(dst))
                setfield!(dst, fname2, getfield(src, fname2))
            end
        end

        backup!(rep)
        setSimCalls!(rep, rep.calls) # requires rep.backup to be defined

        rep.isRunnable = rep.backup.isRunnable = true
    end
end

function simulateRep!(rep::Simulation)
    @assert(rep.isRunnable)
    simulate!(rep)
end

# run simulation replications
# requires replications to be 'runnable', see function makeRepsRunnable!()
# if parallel=true, will use multi-threading
# mutates: reps
function simulateReps!(reps::Vector{Simulation}; parallel::Bool=false, numThreads::Int=maxNumThreads, doPrint::Bool=false)
    @assert(allunique(reps))
    @assert(all(rep -> rep.isRunnable, reps))
    if parallel
        doPrint = false # IO not thread safe
        runParallel!(simulateRep!, reps...; numThreads=numThreads)
    else
        numReps = length(reps)
        for (i, rep) in enumerate(reps)
            doPrint && print("\rSimulating replication $i of $numReps.")
            simulateRep!(rep)
        end
        doPrint && println()
    end
end

function resetRep!(rep::Simulation)
    @assert(rep.isRunnable)
    reset!(rep)
end

function resetReps!(reps::Vector{Simulation}; parallel::Bool=false, numThreads::Int=maxNumThreads)
    @assert(allunique(reps))
    @assert(all(rep -> rep.isRunnable, reps))
    if parallel
        runParallel!(resetRep!, reps...; numThreads=numThreads)
    else
        for rep in reps
            resetRep!(rep)
        end
    end
end
