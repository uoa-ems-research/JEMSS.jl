using Pkg

# Install packages which JEMSS needs to know the version numbers of,
# can only seem to get version numbers for installed packages, using `Pkg.installed()`
Pkg.add("ArchGDAL")
Pkg.add("JuMP")
Pkg.add("LightGraphs")

include("../data/unzip_data.jl")
