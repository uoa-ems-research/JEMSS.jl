using ZipFile

function pathrecursion(path::String, f::Function=x -> x)
    if isfile(path)
        return f(path)
    elseif isdir(path)
        return vcat([pathrecursion(joinpath(path, item), f) for item in readdir(path)]...)
    end
end

function unzip(path::String)
    if isfile(path) && endswith(path, ".zip")
        try
            println("Unzipping: ", path)
            reader = ZipFile.Reader(path)
            for file in reader.files
                open(joinpath(dirname(path), file.name), "w") do f
                    write(f, read(file))
                end
            end
            close(reader)
        catch
            println("Failed to unzip.")
        end
    end
end

cd(@__DIR__) do
    # unpack all zip files
    for folder in ("auckland", "edmonton", "manhattan", "utrecht")
        path = joinpath(pwd(), "cities", folder)
        pathrecursion(path, unzip)
    end

    # filenames = pathrecursion(pwd())
    # dirnames = unique(dirname.(filenames)) # or unique(pathrecursion(pwd(), dirname))
end
