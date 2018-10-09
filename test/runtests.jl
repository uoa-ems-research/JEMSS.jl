using JEMSS
using Base.Test

cd(@__DIR__) do
	isdir("temp") || mkpath("temp")
	include("test_demand.jl")
	include("test_demand_coverage.jl")
end
