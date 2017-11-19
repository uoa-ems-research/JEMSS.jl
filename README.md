# AmbulanceSim.jl

[![Build Status](https://travis-ci.org/samridler/AmbulanceSim.jl.svg?branch=master)](https://travis-ci.org/samridler/AmbulanceSim.jl)

[![Coverage Status](https://coveralls.io/repos/samridler/AmbulanceSim.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/samridler/AmbulanceSim.jl?branch=master)

[![codecov.io](http://codecov.io/github/samridler/AmbulanceSim.jl/coverage.svg?branch=master)](http://codecov.io/github/samridler/AmbulanceSim.jl?branch=master)

## Warning
This package is currently being developed.
Expect bugs, and backward compatibility issues between commits.

## Installation
This package (and another that it requires) is unregistered so you will need to `Pkg.clone` it as follows:
```julia
Pkg.clone("https://github.com/yeesian/ArchGDAL.jl.git")
Pkg.clone("https://github.com/samridler/AmbulanceSim.jl.git")
```

## Usage

### Simulation
```julia
using AmbulanceSim

filename = selectXmlFile()
sim = initSimulation(filename);

simulate!(sim) # run the simulation
printSimStats(sim) # some basic statistics
```
For a Windows OS, the function `selectXmlFile()` will open a file selection dialog using powershell, for other operating systems a prompt will be given to enter the filename string.
The simulation is initialised from an xml configuration file which contains a list of input files (files for ambulances, calls, stations, etc.), and other parameters.

### Animation
```julia
using AmbulanceSim
animate(port = 8001)
```
This will open a local server with port number 8001, though other port numbers may be specified. Open a web browser and enter the URL `localhost:8001`. Firefox and Chrome work, Edge does not work well for this, other browsers have not been tested. Multiple browser windows may use this port. The connection may take a few seconds to be established. You will then be prompted for a simulation configuration filename from Julia (not the browser window).
Once the simulation has been initialised, the browser window, using [Mapbox](https://www.mapbox.com/), should show a map of the region of interest, the ambulances, hospitals, stations, and road arcs. Controlling the animation is done with buttons and text input in a box at the bottom right of the window. The 'show arcs' check-box is provided so that the road arcs can be hidden, which is useful for reducing computation of drawing while the animation is running.

![animation_frame](https://i.imgur.com/qZuS5NZ.png)

Notes on animation:

 - for the timing control, 'speed' is the ratio of animation time to real time; speed of 1 gives real-time (real slow) animation. Calculating the ratio correctly currently requires the input files to have time units in days.
 - input files should have (latitude, longitude) coordinates, these correspond with the y and x fields in the input files

### Examples

The example folder contains example.jl, run this using `include("example.jl")` (run from the containing directory). This will generate some artificial data and save it to files, which can then be used as input to simulate with (example.jl script shows some cases of this) or animate with.

## License
This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.
