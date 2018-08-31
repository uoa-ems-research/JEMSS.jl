# read geographic data from file and apply f to data
function readGeoFile(f::Function, filename::String)
	return ArchGDAL.registerdrivers() do
		ArchGDAL.read(filename) do dataset
			f(dataset)
		end
	end
end

# read raster file using ArchGDAL package, return as custom Raster type
function readRasterFile(rasterFilename::String)
	
	# open raster, get data, close raster
	# only gets first raster band, may change this later
	(geoTransform, z) = readGeoFile(rasterFilename) do dataset
		rasterband = ArchGDAL.getband(dataset, 1)
		return (ArchGDAL.getgeotransform(dataset), ArchGDAL.read(rasterband))
	end
	
	# shorthand:
	(nx, ny) = size(z) # number of cells in x and y directions
	(x1, dxdi, dxdj, y1, dydi, dydj) = geoTransform # dxdi is change in x per change in index i, for z[i,j]; likewise for dxdj, dydi, dydj
	dx = dxdi # width of cells in x direction
	dy = dydj # height of cells in y direction (may be negative)
	
	# data checks
	assert(dxdj == 0 && dydi == 0) # otherwise raster is sloping, so changing x index changes y value, and vice-versa
	assert(dx > 0)
	assert(dy != 0)
	
	# convert data for easier use
	# find x and y vectors to represent raster cell centres, make sure values are increasing
	xMin = x1 + 0.5*dx # we know dx > 0
	# (xMin, dx, z) = (dx > 0) ? (x1 + 0.5*dx, dx, z) : (x1 + (nx-0.5)*dx, -dx, flipdim(z,1))
	(yMin, dy, z) = (dy > 0) ? (y1 + 0.5*dy, dy, z) : (y1 + (ny-0.5)*dy, -dy, flipdim(z,2))
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
		for i = 1:length(v)-1
			assert(isapprox(v[i] + dv, v[i+1]; rtol = eps(Float)))
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
function rasterCellRandLocation(raster::Raster, i::Int, j::Int;
	rng::AbstractRNG = Base.GLOBAL_RNG)
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
	assert(1 <= i && i <= length(raster.z)) # ind2sub does not check this
	return ind2sub(raster.z, i)
end

# generate n random locations for a given raster sampler
function rasterRandLocations(rasterSampler::RasterSampler, n::Int)
	randLocations = Vector{Location}(n)
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
