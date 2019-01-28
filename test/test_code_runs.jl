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

# tests to check that code runs without crashing, but not checking that output is expected

# check that different sim configs can be opened, and sim can run
@testset "sim configs" begin
	@assert(isdir("data/regions/small/1"))
	@info("Testing sim configs")
	simConfigFolder = "data/regions/small/1/sim_configs"
	for configFilename in readdir(simConfigFolder)
		filename = joinpath(pwd(), simConfigFolder, configFilename)
		sim = initSim(filename, doPrint = false);
		simulate!(sim)
		@test true
	end
end

# test the scripts in /example/scripts
if true
	@info("Skipped testing of example scripts") # default is to not run this test set
else
	@testset "example scripts" begin
		scriptsFolder = joinpath(dirname(pathof(JEMSS)), "..", "example", "scripts")
		
		# local search scripts, tested with a very small sim for speed
		cd("data/regions/small/2") do
			runGenConfig("gen_config.xml", overwriteOutputPath = true, doPrint = false)
			isdir("output") || mkdir("output")
			for script in ["deployment_local_search.jl", "nested_comp_table_local_search.jl", "priority_list_local_search.jl"]
				@info(string("Running script: ", script))
				include(joinpath(scriptsFolder, script))
				@test true
			end
			println()
			rm("output", recursive = true)
		end
		
		# deployment_ranking.jl
		cd("data/regions/small/1") do
			isdir("output") || mkdir("output")
			@info("Running script: deployment_ranking.jl")
			include(joinpath(scriptsFolder, "deployment_ranking.jl"))
			rm("output", recursive = true)
			@test true
		end
	end
end
