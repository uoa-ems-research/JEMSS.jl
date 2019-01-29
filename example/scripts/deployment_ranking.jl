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

# Script to generate many different deployments, and rank them by batch mean response times.

using JEMSS
using Random
# using Statistics

# parameters:
const configFilename = "sim_config.xml"
const outputFolder = "output"
const batchMeanResponseTimesFilename = joinpath(outputFolder, "deployments_batch_mean_response_times.csv")
const batchTime = 1.0 # batch by day
const warmUpTime = 1.0 # one day warm-up period

# additional parameters to use if config does not include deployments:
const deploymentsOutFilename = joinpath(outputFolder, "deployments.csv")
numDeployments = 2 # number of deployments to generate and test
deploymentRng = Random.MersenneTwister(0)

# load sim
sim = initSim(configFilename, createBackup = true, doPrint = false)
(numAmbs, numStations) = (sim.numAmbs, sim.numStations)
@assert(numAmbs > 0 && numStations > 1) # otherwise only one deployment is possible

# read deployments file
if haskey(sim.inputFiles, "deployments")
	(deployments, numStations2) = readDeploymentsFile(sim.inputFiles["deployments"])
	numDeployments = length(deployments)
	@assert(numStations == numStations2)
elseif haskey(sim.inputFiles, "deploymentPolicies") # compat
	(deployments, numStations2) = readDeploymentPoliciesFile(sim.inputFiles["deploymentPolicies"])
	numDeployments = length(deployments)
	@assert(numStations == numStations2)
else
	# generate unique random deployments
	deployments = makeRandDeployments(numAmbs, numStations, numDeployments; rng = deploymentRng)
	writeDeploymentsFile(deploymentsOutFilename, deployments, numStations)
	println("'deployments' input file not found in config, generated $numDeployments deployments and saved to $deploymentsOutFilename")
end

# run sim for one day, to compile simulation functions
applyDeployment!(sim, deployments[1])
simulateToTime!(sim, 1.0)
resetSim!(sim)

# run simulation for each deployment, save call response times
getCallResponseTimes(sim::Simulation) = [call.responseTime * 24 * 60 for call in sim.calls] # times in minutes
responseTimeUnits = "minutes"
t = time()
callResponseTimes = simulateDeployments!(sim, deployments, getCallResponseTimes; showEta = true)
println("total time for simulating $numDeployments deployments (seconds): ", round(time() - t, digits = 2))

# calculate batch mean response times
x = Vector{Vector{Float}}() # will populate with batch mean response times
callArrivalTimes = [call.arrivalTime for call in sim.calls] # same for each simulation
startTime = sim.startTime + warmUpTime; endTime = callArrivalTimes[end] # ignore last call (cool-down period)
for i = 1:numDeployments
	push!(x, calcBatchMeans(callResponseTimes[i], callArrivalTimes, batchTime;
		startTime = startTime, endTime = endTime, rmPartialBatch = true))
end
batchMeanResponseTimes = collect(hcat(x...)') # batchMeanResponseTimes[i,j] is for calls in simulation i, batch j

@warn("Have assumed that batch time of $batchTime is sufficient for values in batchMeanResponseTimes[i,:] to be from a normal distribution (for each i).")

# Check for serial autocorrelation of batchMeanResponseTimes for each deployment,
# if autocorrelation is detected, then the batch sizes / durations need to be increased.
# Will apply AR(0) model and use Durbin-Watson test.
x = batchMeanResponseTimes # shorthand
dwPValues = [calcAR0DurbinWatsonTestPValue(x[i,:]) for i = 1:size(x,1)]
for p in [0.01, 0.05, 0.10]
	println(" number of p-values <= ", p, ": ", count(dwPValues .<= p), " out of ", numDeployments)
end
# number of p-values <= p should be approximately <= p * numDeployments, for there to be no evidence against null hypothesis (H0: no serial autocorrelation)

# save batch mean response times to file
writeBatchMeanResponseTimesFile(batchMeanResponseTimesFilename, batchMeanResponseTimes;
	batchTime = batchTime, startTime = startTime, endTime = endTime, responseTimeUnits = responseTimeUnits)
println("Saved batch mean response times to $batchMeanResponseTimesFilename")

# calculate mean and standard error of batch mean response times, plot for each deployment
false && begin
using Plots
x = batchMeanResponseTimes # shorthand
conf = 0.95
order = sortperm(squeeze(mean(x, dims = 2),2)) # sort by average value
Plots.plotly()
plot = meanErrorPlot(x[order,:], conf);
Plots.title!(plot, "Deployment : batch mean response time, mean & error");
Plots.xaxis!(plot, "Deployment");
Plots.yaxis!(plot, "Batch mean response time ($responseTimeUnits)");
display(plot);
end
