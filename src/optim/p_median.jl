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

# P-median problem

pMedianDefaultOptions = Dict([:x_bin => true, :y_bin => false, :solver => "gurobi", :solver_args => [], :solver_kwargs => []])
try Gurobi catch; pMedianDefaultOptions[:solver] = "cbc" end

"""
	solvePMedian(c::Array{Float,2}, n::Int; options::Dict{Symbol,T} = pMedianOptions)
Solves the p-median problem for costs `c` and `n` facilities. `c[i,j]` is the cost of serving demand point `j` from facility location `i`.
The p-median problem is to locate facilities in order to minimise the total cost of serving demand, where each demand point is served by the facility which can serve the point at lowest cost.

Returns a dictionary with the objective value and decision variable values.
"""
function solvePMedian(c::Array{Float,2}, n::Int; options::Dict{Symbol,T} = pMedianDefaultOptions) where T <: Any
	@assert(n >= 1)
	s, d = size(c) # number of potential facilities, number of demand points
	
	global pMedianDefaultOptions
	defaultOptions = copy(pMedianDefaultOptions)
	merge!(defaultOptions, options)
	
	m = model = Model()
	
	# solver
	solver = options[:solver]
	args = options[:solver_args]
	kwargs = options[:solver_kwargs]
	jump_ge_0_19 = pkgVersions["JuMP"] >= v"0.19"
	if jump_ge_0_19
		if solver == "cbc" set_optimizer(model, with_optimizer(Cbc.Optimizer, logLevel=0, args...; kwargs...))
		elseif solver == "glpk" set_optimizer(model, with_optimizer(GLPK.Optimizer, args...; kwargs...))
		elseif solver == "gurobi" @stdout_silent(set_optimizer(model, with_optimizer(Gurobi.Optimizer, OutputFlag=0, args...; kwargs...)))
		end
	else
		if solver == "cbc" setsolver(model, CbcSolver(args...; kwargs...))
		elseif solver == "glpk" setsolver(model, GLPKSolverMIP(args...; kwargs...))
		elseif solver == "gurobi" @stdout_silent(setsolver(model, GurobiSolver(OutputFlag=0, args...; kwargs...)))
		end
	end
	
	# variables
	options[:x_bin] ? @variable(m, x[i=1:s], Bin) : @variable(m, 0 <= x[i=1:s] <= 1) # x[i] = 1 if one facility at location i, 0 otherwise
	options[:y_bin] ? @variable(m, y[i=1:s,j=1:d], Bin) : @variable(m, 0 <= y[i=1:s,j=1:d] <= 1) # y[i,j] = 1 if facility at location i serves point j
	
	@constraints(model, begin
		(maxFacilityCount, sum(x) <= n)
		(serveDemandFromOpenFacilities[i=1:s,j=1:d], y[i,j] <= x[i]) # can only serve point j from facility location i if it is opened
		(demandServedBySingleFacility[j=1:d], sum(y[:,j]) == 1)
	end)
	
	@expressions(model, begin
		cost, sum(c.*y)
	end)
	
	# solve
	if jump_ge_0_19
		@objective(model, Min, cost)
		@stdout_silent optimize!(model)
		@assert(termination_status(model) == MOI.OPTIMAL)
	else
		@objective(model, :Min, cost)
		status = @stdout_silent solve(model)
		@assert(status == :Optimal)
	end
	
	# get solution
	results = Dict()
	results[:cost] = jump_ge_0_19 ? JuMP.value.(cost) : getvalue(cost)
	results[:x] = jump_ge_0_19 ? JuMP.value.(x) : getvalue(x)
	results[:y] = jump_ge_0_19 ? JuMP.value.(y) : getvalue(y)
	results[:model] = model
	
	return results
end
