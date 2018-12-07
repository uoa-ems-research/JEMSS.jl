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

# crop raster borders to fit in map borders, return true if successful, false otherwise
# mutates: raster
function cropRaster!(raster::Raster, map::Map)
	
	# find x, y, and z values to keep
	ix = find((map.xMin + raster.dx / 2 .<= raster.x) .* (raster.x .<= map.xMax - raster.dx / 2))
	iy = find((map.yMin + raster.dy / 2 .<= raster.y) .* (raster.y .<= map.yMax - raster.dy / 2))
	x = raster.x[ix]
	y = raster.y[iy]
	z = raster.z[ix,iy]
	
	# create cropped raster
	croppedRaster = nothing # init
	try
		croppedRaster = Raster(x, y, z)
	catch
		println("Raster does not overlap enough with map; raster not cropped")
		@show length(x), length(y)
		return false
	end
	
	# overwrite raster field values in raster with those in croppedRaster
	for fname in fieldnames(raster)
		setfield!(raster, fname, getfield(croppedRaster, fname))
	end
	
	return true
end

# given raster and x and y indices (i, j) for a cell,
# generate a random location within the cell
# (x[i] +/- 0.5*dx, y[j] +/- 0.5*dy)
function rasterCellRandLocation(raster::Raster, i::Int, j::Int;
	rng::AbstractRNG = GLOBAL_RNG)
	randLocation = Location()
	randLocation.x = raster.x[i] + (rand(rng) - 0.5) * raster.dx
	randLocation.y = raster.y[j] + (rand(rng) - 0.5) * raster.dy
	return randLocation
end

# function rasterZIndexToXIndex(raster::Raster, i::Int)
	# # return index of raster.x for given index of raster.z
	# return mod(i-1, raster.nx) + 1 # x index
# end

# function rasterZIndexToYIndex(raster::Raster, i::Int)
	# # return index of raster.y for given index of raster.z
	# return cld(i, raster.nx) # y index
# end

function rasterZIndexToXYIndices(raster::Raster, i::Int)
	# return (x,y) indices for raster.z[i]
	@assert(1 <= i && i <= length(raster.z)) # ind2sub does not check this
	return ind2sub(raster.z, i)
end

# generate n random locations for a given raster sampler
function rasterRandLocations(rasterSampler::RasterSampler, n::Int)
	randLocations = Vector{Location}(undef, n)
	zis = rand(rasterSampler.cellDistrRng, n) # z indices of random locations
	for i = 1:n
		zi = zis[i]
		(xi, yi) = rasterZIndexToXYIndices(rasterSampler.raster, zi)
		randLocations[i] = rasterCellRandLocation(rasterSampler.raster, xi, yi; rng = rasterSampler.cellLocRng)
	end
	return randLocations
end

function printRasterSize(raster::Raster)
	println("nx * ny = ", raster.nx * raster.ny)
	println("xMin; xMax; yMin; yMax:")
	println(minimum(raster.x) - raster.dx/2, "; ", maximum(raster.x) + raster.dx/2, "; ",
		minimum(raster.y) - raster.dy/2, "; ", maximum(raster.y) + raster.dy/2)
end
