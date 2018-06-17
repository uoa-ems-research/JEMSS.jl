# JEMSS.jl

[![Build Status](https://travis-ci.org/samridler/JEMSS.jl.svg?branch=master)](https://travis-ci.org/samridler/JEMSS.jl)

[![Coverage Status](https://coveralls.io/repos/samridler/JEMSS.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/samridler/JEMSS.jl?branch=master)

[![codecov.io](http://codecov.io/github/samridler/JEMSS.jl/coverage.svg?branch=master)](http://codecov.io/github/samridler/JEMSS.jl?branch=master)

## Warning
This package is currently being developed.
Expect bugs, and backward compatibility issues between commits.

## Installation
This package (and others that it requires) is unregistered so you will need to `Pkg.clone` it as follows:
```julia
Pkg.clone("https://github.com/visr/GDAL.jl.git")
Pkg.clone("https://github.com/yeesian/ArchGDAL.jl.git")
Pkg.clone("https://github.com/samridler/JEMSS.jl.git")
```

## Usage

### Simulation
```julia
using JEMSS

filename = selectXmlFile()
sim = initSimulation(filename);

simulate!(sim) # run the simulation
printSimStats(sim) # some basic statistics
```
For a Windows OS, the function `selectXmlFile()` will open a file selection dialog using powershell, for other operating systems a prompt will be given to enter the filename string.
The simulation is initialised from an xml configuration file which contains a list of input files (files for ambulances, calls, stations, etc.), and other parameters.

### Animation
```julia
using JEMSS
animate(port = 8001, configFilename = filename)
```
This will open a web browser window to `localhost:8001` (other port numbers may be used) and load the simulation for the given configuration file.
The connection may take a few seconds to be established.
Once the simulation has been initialised, the browser window, using [Mapbox](https://www.mapbox.com/), will show a map of the region of interest, the ambulances, hospitals, stations, and road arcs.
Controlling the animation is done with buttons and text input in a box at the bottom right of the window.
The 'show arcs' check-box is provided so that the road arcs can be hidden, which is useful for reducing computation of drawing while the animation is running.

![animation_frame](https://i.imgur.com/CN6Yeak.png)

Notes on animation:

- Firefox and Chrome work, Edge does not work well, other browsers have not been tested.
- Multiple browser windows may use the same port.
- If you do not specify a simulation configuration filename, then you will be prompted for this from Julia (not the browser window).
- For the timing control, 'speed' is the ratio of animation time to real time; speed of 1 gives real-time (real slow) animation. Calculating the ratio correctly currently requires the input files to have time units in days.
- Input files should have (latitude, longitude) coordinates, these correspond with the y and x fields in the input files.

### Examples

The example folder contains example.jl, run this using `include("example.jl")` (run from the containing directory). This will generate some artificial data and save it to files, which can then be used as input to simulate with (example.jl script shows some cases of this) or animate with.

## License
This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.
