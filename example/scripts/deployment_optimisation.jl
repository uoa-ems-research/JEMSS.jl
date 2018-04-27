# Script to generate many different deployment policies, and rank them by batch mean response times.
# Run from julia:
# ARGS = [simConfigFilename, numDepols] # these are optional; depol = deployment policy
# include("deployment_optimisation.jl")

using JEMSS

# load sim
configFilename = (length(ARGS) >= 1 ? ARGS[1] : selectXmlFile())
sim = initSimulation(configFilename)
(numAmbs, numStations) = (length(sim.ambulances), length(sim.stations))
assert(numAmbs > 0 && numStations > 1) # otherwise only one deployment policy is possible

# read deployment policies file
if haskey(sim.inputFiles, "deploymentPolicies")
	(depols, numStations2) = readDeploymentPoliciesFile(sim.inputFiles["deploymentPolicies"].path)
	numDepols = length(depols)
	assert(numStations == numStations2)
else
	# generate unique random deployment policies
	numDepols = (length(ARGS) >= 2 ? ARGS[2] : 100) # number of deployment policies to test
	depols = makeRandDeploymentPolicies(numAmbs, numStations, numDepols)
	filename = joinpath(sim.outputPath, "deployment_policies.csv")
	writeDeploymentPoliciesFile(filename, depols, numStations)
	println("'deploymentPolicies' input file not found in config, generated $numDepols deployment policies and saved to $filename")
end

# run sim for one day, to compile simulation functions
applyDeploymentPolicy!(sim, depols[1])
simulateToTime!(sim, 1.0)
resetSim!(sim)

# run simulation for each deployment policy, save call response times
getCallResponseTimes(sim::Simulation) = [call.responseTime * 24 * 60 for call in sim.calls] # times in minutes
responseTimeUnits = "minutes"
t = time()
callResponseTimes = simulateDeploymentPolicies!(sim, depols, getCallResponseTimes; showEta = true)
println("total time for simulating $numDepols deployment policies (seconds): ", round(time() - t, 2))

# calculate batch mean response times
x = [] # will populate with batch mean response times
callArrivalTimes = [call.arrivalTime for call in sim.calls] # same for each simulation
batchTime = 7.0; startTime = sim.startTime + 1.0; endTime = callArrivalTimes[end] # batch response times by week, with one day warm-up; ignore last call (cool-down period)
for i = 1:numDepols
	push!(x, calcBatchMeans(values = callResponseTimes[i], times = callArrivalTimes,
		batchTime = batchTime, startTime = startTime, endTime = endTime)[1])
end
batchMeanResponseTimes = hcat(x...)' # batchMeanResponseTimes[i,j] is for calls in simulation i, batch j

warn("Have assumed that batch time of $batchTime is sufficient for values in batchMeanResponseTimes[i,:] to be from a normal distribution (for each i).")

# Check for serial autocorrelation of batchMeanResponseTimes for each policy,
# if autocorrelation is detected, then the batch sizes / durations need to be increased.
# Will apply AR(0) model and use Durbin-Watson test.
x = batchMeanResponseTimes # shorthand
dwPValues = [calcAR0DurbinWatsonTestPValue(x[i,:]) for i = 1:size(x,1)]
for p in [0.01, 0.05, 0.10]
	println(" number of p-values <= ", p, ": ", count(dwPValues .<= p), " out of ", numDepols)
end
# number of p-values <= p should be approximately <= p * numDepols, for there to be no evidence against null hypothesis (H0: no serial autocorrelation)

# save batch mean response times to file
filename = joinpath(sim.outputPath, "depols_batch_mean_response_times.csv")
writeBatchMeanResponseTimesFile(filename, batchMeanResponseTimes;
	batchTime = batchTime, startTime = startTime, endTime = endTime, responseTimeUnits = responseTimeUnits)
println("Saved batch mean response times to $filename")

# calculate mean and standard error of batch mean response times, plot for each deployment policy
x = batchMeanResponseTimes # shorthand
conf = 0.95
order = sortperm(squeeze(mean(x,2),2)) # sort by average value
Plots.plotly()
plot = meanErrorPlot(x[order,:], conf);
Plots.title!(plot, "Deployment policy : batch mean response time, mean & error");
Plots.xaxis!(plot, "Deployment policy");
Plots.yaxis!(plot, "Batch mean response time ($responseTimeUnits)");
display(plot);
