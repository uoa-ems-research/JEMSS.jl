# JEMSS.jl
Julia package for Emergency Medical Services Simulation

<!-- [![Build Status](https://travis-ci.com/uoa-ems-research/JEMSS.jl.svg?branch=master)](https://travis-ci.com/uoa-ems-research/JEMSS.jl) -->
<!-- [![Coverage Status](https://coveralls.io/repos/github/uoa-ems-research/JEMSS.jl/badge.svg?branch=master)](https://coveralls.io/github/uoa-ems-research/JEMSS.jl?branch=master) -->
<!-- [![codecov.io](http://codecov.io/github/uoa-ems-research/JEMSS.jl/coverage.svg?branch=master)](http://codecov.io/github/uoa-ems-research/JEMSS.jl?branch=master) -->

## Installation
To install and build this package, run the following commands in the Pkg REPL mode (entered by pressing `]` from the Julia REPL; press backspace to get back):
```
pkg> add https://github.com/uoa-ems-research/JEMSS.jl
pkg> build JEMSS
```
The build step unpacks zipped files in the data folder.

## Simulation example
To run an example script that loads and runs a simulation and then writes statistics to files:
```julia
using JEMSS
include(joinpath(JEMSS.jemssDir, "example/example.jl"))
```
This [example](example/example.jl) is for Auckland city, focusing on the urban area.
The simulation is initialised from a configuration file ([sim_config.xml](example/input/sim_config.xml)) which contains a list of input files (files for ambulances, calls, stations, road network, etc.) and other parameters.
The output files will be written to the folder example/output.
The first time that this script is run it will take an additional minute or so to compute and serialise the all-pairs shortest-path data for the road network; subsequent runs will be faster as they read the serialised data.
After this script has been run, the simulation can be animated with `animate!(sim)` provided that the [animation setup](#animation-setup) steps have been completed.

A list of further example scripts that may be useful can be found in [example/other_examples.jl](example/other_examples.jl).

## Cities
The package includes these city (/island/region) models: Auckland, New Zealand; Edmonton, Alberta, Canada; Manhattan Island, New York, USA; and Utrecht, Netherlands.
Example scripts to simulate these cities can be found in [example/other_examples.jl](example/other_examples.jl).
The models are not exact replicas of the cities, but they are at least city-like.
More information on the city models (sources, simulation results) is available in this [pdf document (on Google drive)](https://tinyurl.com/2p8m9a8p).

## Animation

### Setup
In order to use the animation tool, you will need a mapbox access token.
First create an account at [mapbox.com](https://www.mapbox.com/) then create an access token, the default settings should be sufficient.
Then copy the token and:
- For JEMSS version >= 1.3.6:
Call `JEMSS.setMapboxAccessToken("your token here")` which will save the token to JEMSS/src/animation/mapbox_access_token.txt.
- For JEMSS version < v1.3.6:
Paste the token into JEMSS/src/animation/index.html, search for and replace the value of `L.mapbox.accessToken`.

You should now be able to use the animation tool.
As of writing (June 2024) the use of mapbox is free up to 200,000 tile requests per month which should be plenty; I average under 1,000 per month.

I had previously put my own mapbox access token in the open source code but this is a big no-no.
I have now deleted this token from my mapbox account, if the animation tool no longer shows the map then this may be the reason, my apologies.

### Run
To animate a simulation:
```julia
using JEMSS
sim = initSim("config_filename");
animate!(sim; port = 8001)
```
The call to `animate!` will open a web browser window to `localhost:8001` (other port numbers may be used).
The connection may take a few seconds to be established.
The browser window, using [Mapbox](https://www.mapbox.com/), will show a map containing the `sim` region, the ambulances, hospitals, stations, and roads.
Controlling the animation is done with buttons and text input in a box at the bottom right of the window.
To have lines drawn between each ambulance and its destination, check the 'Show destinations' box.
The 'Show road arcs' check-box is provided so that the road network arcs can be hidden, which reduces the computation required to display the city while the simulation is running.

![animation_frame](assets/animation/animation-frame.png)

Notes on animation:

- Firefox and Chrome work, Edge does not work well, other browsers have not been tested.
- Multiple browser windows may use the same port.
- For the timing control, 'Sim speed' is the ratio of simulation time to real time; speed of 1 gives real-time (real slow) simulation. This requires the input files to have time units in days.
- Input files should have (latitude, longitude) coordinates, these correspond with the y and x fields in the input files.

## Misc
For solving linear and integer programs, CBC and GLPK solvers are used, though [Gurobi](https://www.gurobi.com/) will be used (for difficult problems such as p-median and DDSM in JEMSS) if it is installed along with the [Gurobi.jl](https://github.com/jump-dev/Gurobi.jl) package, as Gurobi generally solves faster.

Backslashes are special characters in Julia strings and so if a path includes backslashes (e.g. `"path\to\file.txt"`), it needs to be handled as a raw string (`raw"path\to\file.txt"`).

### Adding cities
If you have a city model that you would like to add, please submit a pull request.
Cities should follow the same folder structure as those existing (in [data/cities](data/cities)), with a data folder containing the raw data along with any sources and licenses, and a model folder containing the input files for the calls, stations, road network, etc.
Note that any large files should be compressed as a .zip file.

## License
This project is licensed under the Apache License 2.0 - see the [LICENSE.md](LICENSE.md) file for details.
