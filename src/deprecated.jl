
## renamed "deployment policy" to "deployment"

const Depol = Deployment # not sure how to deprecate a const

@deprecate(
makeRandDeploymentPolicy(numAmbs::Int, numStations::Int; rng::AbstractRNG = Base.GLOBAL_RNG),
makeRandDeployment(numAmbs, numStations; rng = rng))

@deprecate(
makeRandDeploymentPolicies(numAmbs::Int, numStations::Int, numDeployments::Int; rng::AbstractRNG = Base.GLOBAL_RNG),
makeRandDeployments(numAmbs, numStations, numDeployments; rng = rng))

@deprecate(
applyDeploymentPolicy!(sim::Simulation, deployment::Deployment),
applyDeployment!(sim, deployment))

@deprecate(
simulateDeploymentPolicy!(sim::Simulation, deployment::Deployment),
simulateDeployment!(sim, deployment))

@deprecate(
simulateDeploymentPolicies!(sim::Simulation, deployments::Vector{Deployment}, f::Function; showEta::Bool = false),
simulateDeployments!(sim, deployments, f; showEta = showEta))

# write deployments to file
function writeDeploymentPoliciesFile(filename::String, deployments::Vector{Deployment}, numStations::Int)
	warn("'writeDeploymentPoliciesFile' is deprecated, should use 'writeDeploymentsFile' instead")
	numAmbs = length(deployments[1])
	@assert(numStations >= maximum([maximum(d) for d in deployments]))
	numDeployments = length(deployments)
	
	miscTable = Table("miscData", ["numStations", "numDepols"]; rows = [[numStations, numDeployments]])
	deploymentPoliciesTable = Table("deploymentPolicies",
		["ambIndex", ["policy_$i stationIndex" for i = 1:numDeployments]...];
		cols = [collect(1:numAmbs), deployments...])
	writeTablesToFile(filename, [miscTable, deploymentPoliciesTable])
end

function readDeploymentPoliciesFile(filename::String)
	warn("'readDeploymentPoliciesFile' is deprecated, should save deployment with 'writeDeploymentsFile' and load again with 'readDeploymentsFile'")
	tables = readTablesFromFile(filename)
	
	# get counts
	table = tables["miscData"]
	numStations = table.columns["numStations"][1]
	numDeployments = table.columns["numDepols"][1]
	
	# deployments
	table = tables["deploymentPolicies"]
	columns = table.columns # shorthand
	# check number of ambulances
	ambIndexCol = columns["ambIndex"]
	@assert(ambIndexCol == collect(1:length(ambIndexCol)))
	# check that table header has all deployments
	for i = 1:numDeployments
		@assert(in("policy_$i stationIndex", table.header), "Missing deployment $i")
	end
	# create deployments from table
	deployments = Vector{Deployment}(numDeployments)
	for i = 1:numDeployments
		deployments[i] = columns["policy_$i stationIndex"]
		
		for value in deployments[i]
			@assert(in(value, 1:numStations))
		end
	end
	
	return deployments, numStations
end

##
