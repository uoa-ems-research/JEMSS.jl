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

using JEMSS
using Test

cd(@__DIR__) do
    isdir("temp") || mkpath("temp")
    runGenConfig("data/cities/small/1/gen_config.xml", overwriteOutputPath=true, doPrint=false)
    include("test_grid.jl")
    include("test_travel.jl")
    include("test_network.jl")
    include("test_route.jl")
    include("test_demand.jl")
    include("test_demand_coverage.jl")
    include("test_zhang_ip.jl")
    include("test_mclp.jl")
    include("test_nested_mclp.jl")
    include("test_binary_search.jl")
    include("test_histogram.jl")
    include("test_code_runs.jl")
end
