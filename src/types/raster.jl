# read raster file using RasterIO package, return as custom Raster type
# note that RasterIO does not automatically check if raster cells are evenly spaced
function readRasterFile(rasterFilename::String)
	
	# open raster, get data, close raster
	assert(isfile(rasterFilename))
	raster = RasterIO.openraster(rasterFilename, RasterIO.GA_ReadOnly)
	z = RasterIO.fetch(raster,1)
	data = RasterIO.geotransform(raster)
	RasterIO.closeraster(raster)
	
	# shorthand:
	nx = raster.width # number of cells in x direction
	ny = raster.height # number of cells in y direction
	dx = data[2] # width of cells in x direction
	dy = data[6] # height of cells in y direction (may be negative)
	
	# data checks
	assert(raster.nband == 1) # otherwise would have to choose which band to use
	assert(nx == size(z,1) && ny == size(z,2))
	assert(data[3] == 0 && data[5] == 0) # otherwise raster is sloping, so changing x index changes y value, and vice-versa
	assert(dx > 0) # otherwise RasterIO.openraster should have failed
	assert(dy != 0)
	
	# convert data for easier use
	# find x and y vectors to represent raster cell centres, make sure values are increasing
	xMin = data[1] + 0.5*dx # we know dx > 0
	# (xMin, dx, z) = (dx > 0) ? (data[1] + 0.5*dx, dx, z) : (data[1] + (nx-0.5)*dx, -dx, flipdim(z,1))
	(yMin, dy, z) = (dy > 0) ? (data[4] + 0.5*dy, dy, z) : (data[4] + (ny-0.5)*dy, -dy, flipdim(z,2))
	x = collect(range(xMin, dx, nx))
	y = collect(range(yMin, dy, ny))
	
	return Raster(x, y, z)
end

# For reading raster stored in jld file, with x, y, and z variables,
# where x and y values are coordinates (each sorted by increasing value)
# and x[i], y[j] corresponds with cell value z[i,j] ((length(x), length(y)) == size(z))
function readJldRasterFile(rasterFilename::String)
	# open and read file
	assert(isfile(rasterFilename))
	data = JLD.load(rasterFilename)
	raster = Raster(data["x"], data["y"], data["z"])
	
	# check data
	
	# z[i,j] should correspond with x[i], y[j]
	assert((length(raster.x), length(raster.y)) == size(raster.z))
	
	# x and y values should increase with index
	assert(issorted(raster.x, lt = <))
	assert(issorted(raster.y, lt = <))
	
	# check that x and y values have (roughly) constant step sizes
	for (v, dv) in [(raster.x, raster.dx), (raster.y, raster.dy)]
		if length(v) > 1
			assert(maximum(abs(v[1:end-1] - (v[2:end] - dv))) .<= eps() * maximum(abs(v)))
		end
	end
	
	return raster
end

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
function rasterCellRandLocation(raster::Raster, i::Int, j::Int)
	randLocation = Location()
	randLocation.x = raster.x[i] + (rand() - 0.5) * raster.dx
	randLocation.y = raster.y[j] + (rand() - 0.5) * raster.dy
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
	assert(1 <= i && i <= length(raster.z)) # ind2sub does not check this
	return ind2sub(raster.z, i)
end

# given a raster where z values are proportional to location (x,y) probabilities
# for each raster cell, return n (=numLocations) many random locations,
# each uniformly spatially distributed within the corresponding cell
function rasterRandLocations(raster::Raster, numLocations::Int)
	
	# check that probabilities are non-negative
	for z in raster.z
		assert(z >= 0)
	end
	
	# cumulative distribution function for z
	zcdf = cumsum(raster.z[:])
	assert(zcdf[end] > 0)
	zcdf /= zcdf[end] # scale to 1
	# zcdf[i] = probability of generating a location within cells for z[1:i]
	
	# # convert x, y data to format where z[i,j] is for x[i,j], y[i,j]
	# x = repmat(raster.x, 1, raster.ny)
	# y = repmat(raster.y', raster.nx, 1)
	
	# use inverse transform sampling to randomly select cells
	u = sort(1-rand(numLocations)) # sorted uniform samples from (0,1], to be matched to zcdf
	randOrder = randperm(numLocations) # u[i] is for location[randOrder[i]]
	randLocations = Vector{Location}(numLocations)
	zi = 1 # current index in z
	for i = 1:numLocations
		# set location of randLocations[randOrder[i]]
		while u[i] > zcdf[zi] || raster.z[zi] == 0
			zi += 1
		end
		# put location somewhere in cell zi
		(xi, yi) = rasterZIndexToXYIndices(raster, zi)
		randLocations[randOrder[i]] = rasterCellRandLocation(raster, xi, yi) # do not use randLocations[i] = location
	end
	
	return randLocations
end

function printRasterSize(raster::Raster)
	println("nx * ny = ", raster.nx * raster.ny)
	println("xMin; xMax; yMin; yMax:")
	println(minimum(raster.x) - raster.dx/2, "; ", maximum(raster.x) + raster.dx/2, "; ",
		minimum(raster.y) - raster.dy/2, "; ", maximum(raster.y) + raster.dy/2)
end
