# JEMSS.jl

[![Build Status](https://travis-ci.com/uoa-ems-research/JEMSS.jl.svg?branch=master)](https://travis-ci.com/uoa-ems-research/JEMSS.jl)
[![Coverage Status](https://coveralls.io/repos/github/uoa-ems-research/JEMSS.jl/badge.svg?branch=master)](https://coveralls.io/github/uoa-ems-research/JEMSS.jl?branch=master)
[![codecov.io](http://codecov.io/github/uoa-ems-research/JEMSS.jl/coverage.svg?branch=master)](http://codecov.io/github/uoa-ems-research/JEMSS.jl?branch=master)

## Warning
This package is currently being developed.
Expect bugs, and backward compatibility issues between commits.

## Installation
To install this package, run the following command in the Pkg REPL mode:
```julia
pkg> add https://github.com/uoa-ems-research/JEMSS.jl
```

## Usage

### Simulation
```julia
using JEMSS

filename = selectXmlFile()
sim = initSim(filename);

simulate!(sim) # run the simulation
printSimStats(sim) # some basic statistics
```
For a Windows OS, the function `selectXmlFile()` will open a file selection dialog using powershell, for other operating systems a prompt will be given to enter the filename string.
The simulation is initialised from an xml configuration file which contains a list of input files (files for ambulances, calls, stations, etc.), and other parameters.

### Animation
```julia
using JEMSS
animate(configFilename = filename, port = 8001)
# alternative:
sim = initSim(filename);
animate!(sim) # will use same port as before (8001)
```
The call to `animate` will open a web browser window to `localhost:8001` (other port numbers may be used) and load the simulation for the given configuration file.
If a simulation is already loaded beforehand then `animate!(sim)` can be used instead.
The connection may take a few seconds to be established.
Once the simulation has been initialised, the browser window, using [Mapbox](https://www.mapbox.com/), will show a map of the region of interest, the ambulances, hospitals, stations, and road arcs.
Controlling the animation is done with buttons and text input in a box at the bottom right of the window.
To have lines drawn between each ambulance and its destination, check the 'Show routes' box.
The 'Show arcs' check-box is provided so that the road arcs can be hidden, which is useful for reducing computation of drawing while the animation is running.

![animation_frame](https://i.imgur.com/UsUnPIw.png)

Notes on animation:

- Firefox and Chrome work, Edge does not work well, other browsers have not been tested.
- Multiple browser windows may use the same port.
- If you do not specify a simulation or configuration filename, then you will be prompted for a configuration filename from Julia (not the browser window).
- For the timing control, 'Sim speed' is the ratio of simulation time to real time; speed of 1 gives real-time (real slow) simulation. This requires the input files to have time units in days.
- Input files should have (latitude, longitude) coordinates, these correspond with the y and x fields in the input files.

### Examples

The example folder contains example.jl, run this using `include("example.jl")` (run from the containing directory). This will generate some artificial data and save it to files, which can then be used as input to simulate with (example.jl script shows some cases of this) or animate with.

## License
This project is licensed under the Apache License 2.0 - see the [LICENSE.md](LICENSE.md) file for details.
