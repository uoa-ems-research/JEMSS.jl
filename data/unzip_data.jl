using BinaryProvider

# Zip files should have been created using BinaryProvider.package(), otherwise unzipping may result in an error.

function pathrecursion(path::String, f::Function = x->x)
	if isfile(path)
		return f(path)
	elseif isdir(path)
		return vcat([pathrecursion(joinpath(path, item), f) for item in readdir(path)]...)
	end
end

function unzip(path::String)
	if isfile(path) && endswith(path, ".7z")
		println("Unzipping: ", path)
		unpack(path, dirname(path))
	end
end

cd(@__DIR__) do
	# unpack all zip files
	for folder in ("auckland",)
		path = joinpath(pwd(), "cities", folder)
		pathrecursion(path, unzip)
	end
	
	# filenames = pathrecursion(pwd())
	# dirnames = unique(dirname.(filenames)) # or unique(pathrecursion(pwd(), dirname))
end
