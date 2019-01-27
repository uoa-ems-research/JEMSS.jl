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
	simConfigFolder = "data/regions/small/1/sim_configs"
	for configFilename in readdir(simConfigFolder)
		@info(string("Running sim config: ", configFilename))
		filename = joinpath(pwd(), simConfigFolder, configFilename)
		sim = initSim(filename, doPrint = false);
		simulate!(sim)
		@test true
	end
end
