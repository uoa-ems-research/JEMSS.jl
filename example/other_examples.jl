@warn("This script is not intended to be run, it is meant as a reference for further example scripts.")

# Auckland, with single replication.
# This is the same as the example script: JEMSS/example/example.jl.
# To animate (after sim is loaded), run: animate!(sim)
include("$(JEMSS.jemssDir)/example/cities/auckland/urban/1/single/example.jl")

# Auckland, with multiple replications.
# Multiple replications cannot be animated.
include("$(JEMSS.jemssDir)/example/cities/auckland/urban/1/reps/example.jl")

# Edmonton, with multiple replications.
include("$(JEMSS.jemssDir)/example/cities/edmonton/1/reps/example.jl")

# Manhattan, with multiple replications.
include("$(JEMSS.jemssDir)/example/cities/manhattan/1/reps/example.jl")

# Small and artificially generated example; sim loads relatively quickly.
include("$(JEMSS.jemssDir)/example/cities/generated/small/1/example.jl")
