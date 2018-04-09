# misc file input/output functions

Base.transpose(s::String) = s # not sure where else to put this

function readDlmFileNextLine!(file::IOStream; delim::Char = delimiter)
	try
		return readdlm(IOBuffer(readline(file)), delim)
	catch
		return []
	end
end

function readDlmFile(filename::String; delim::Char = delimiter)
	@assert(isfile(filename), "file not found: $filename")
	file = open(filename, "r")
	data = readdlm(file, delim)
	close(file)
	return data
end

# open file for writing after checking that it does not already exist
function openNewFile(filename::String)
	@assert(!isfile(filename), "file already exists: $filename")
	file = open(filename, "w")
	return file
end

# write a line to file, add delimiter between args
function writeDlmLine!(file::IOStream, args...; delim::Char = delimiter)
	for arg in args
		write(file, string(arg), delim)
	end
	write(file, newline)
	# writedlm(file, hcat(args...), delim) # several times slower
end

# create dict with array entries as keys, indices as values
function arrayDict(array::Vector{Any})
	# return Dict(array[i] => i for i = 1:length(array)) # slower for some reason?
	dict = Dict()
	for i = 1:length(array)
		dict[array[i]] = i
	end
	return dict
end

immutable Table
	name::String
	header::Vector{Any} # header names for columns in data, header[i] is for data[:,i]
	headerDict::Dict{Any,Int} # headerDict[h] gives index of h in header
	data::Array{Any,2}
	# columns::Dict{Any,SubArray{Any,1}} # columns[h] gives subarray column data for header h
	columns::Dict{Any,Vector{Any}}
	# columns[h] is same as data[:,headerDict[h]]
	
	function Table(name, header, data)
		name = convert(String, name)
		header = convert(Vector{Any}, header)
		data = convert(Array{Any,2}, data)
		
		(numRows, numCols) = size(data)
		assert(length(header) == numCols)
		assert(allunique(header))
		
		headerDict = arrayDict(header)
		columns = Dict{Any,Vector{Any}}()
		for i = 1:numCols
			columns[header[i]] = data[:,i]
			# columns[header[i]] = view(data, :, i)
		end
		
		return new(name, header, headerDict, data, columns)
	end
	
	function Table(name, header; rows = [], cols = [])
		assert(isempty(rows) + isempty(cols) == 1)
		if !isempty(rows)
			return Table(name, header, ctranspose(hcat(rows...)))
		elseif !isempty(cols)
			return Table(name, header, hcat(cols...))
		end
	end
end

function writeTablesToFile!(file::IOStream, table::Table;
	writeNumRows::Bool = true, writeNumCols::Bool = false)
	
	(numRows, numCols) = size(table.data)
	writeDlmLine!(file, table.name, writeNumRows ? numRows : "", writeNumCols ? numCols : "")
	writeDlmLine!(file, table.header...)
	writedlm(file, table.data, delimiter)
	writeDlmLine!(file, "") # newline
end

function writeTablesToFile!(file::IOStream, tables::Vector{Table})
	for table in tables
		writeTablesToFile!(file, table)
	end
end

function writeTablesToFile!(file::IOStream, tables::Dict{String, Table})
	for (name, table) in tables
		writeTablesToFile!(file, table)
	end
end

function writeTablesToFile(filename::String, tables::Vector{Table})
	file = open(filename, "w")
	writeTablesToFile!(file, tables)
	close(file)
end

function writeTablesToFile(filename::String, table::Table)
	writeTablesToFile(filename, [table])
end

function writeTablesToFile(filename::String, tables::Dict{String, Table})
	file = open(filename, "w")
	writeTablesToFile!(file, tables)
	close(file)
end

function readTablesFromFile(filename::String)
	return readTablesFromData(readDlmFile(filename))
end

function readTablesFromData(data::Array{Any,2})
	(numDataRows, numDataCols) = size(data)
	
	tables = Dict{String, Table}()
	i = 1
	while i <= numDataRows
		# find start of table
		while i <= numDataRows && isempty(data[i,1])
			i += 1
		end
		if i > numDataRows
			break
		end
		
		tableName = data[i,1]
		
		# read row and column counts, if given
		isempty(data[i,2]) ? numRows = nullIndex : numRows = data[i,2]
		isempty(data[i,3]) ? numCols = nullIndex : numCols = data[i,3]
		
		i += 1 # move to table header
		
		# determine column count
		j = 1
		while j <= numDataCols && !isempty(data[i,j])
			j += 1
		end
		numCols == nullIndex ? numCols = j - 1 : assert(numCols == j - 1)
		
		# read table header
		tableHeader = data[i,1:numCols]
		
		i += 1 # move to table data
		
		# determine row count
		j = i
		while j <= numDataRows && !isempty(data[j,1])
			j += 1
		end
		numRows == nullIndex ? numRows = j - i : assert(numRows == j - i)
		
		# create table
		tableData = data[i-1+(1:numRows), 1:numCols]
		table = Table(convert(String,tableName), tableHeader, tableData)
		
		# add table to tables
		assert(!haskey(tables, table.name))
		tables[table.name] = table
		
		# go to next table
		i += numRows # i = j
		assert(i > numDataRows || isempty(data[i,1]))
	end
	
	return tables
end

# crc checksum
function fileChecksum(filename::String)
	return Base.crc32c(Mmap.mmap(filename))
end

function serializeToFile(filename::String, data::Any)
	file = open(filename, "w")
	serialize(file, data)
	close(file)
end

function deserializeFile(filename::String)
	file = open(filename, "r")
	data = deserialize(file)
	close(file)
	return data
end

function interpolateString(s::String)
	return eval(parse(string("\"", s, "\"")))
end

# some convenient functions for reading xml files
xmlFileRoot(filename::String) = root(parse_file(filename))
findElt = find_element # shorthand
eltContent(parentElt::XMLElement, eltString::String) = try content(findElt(parentElt, eltString));
	catch error("Element not found: $eltString"); end
eltContentVal(parentElt::XMLElement, eltString::String) = eval(parse(eltContent(parentElt, eltString)))
eltContentInterpVal(parentElt::XMLElement, eltString::String) = interpolateString(eltContent(parentElt, eltString))

function childrenNodeNames(parentElt::XMLElement)
	childNodes = Vector{String}(0)
	for childNode in child_nodes(parentElt)
		if is_elementnode(childNode)
			push!(childNodes, name(childNode))
		end
	end
	return childNodes
end

function selectXmlFile(; message::String = "Enter xml filename: ")
	if is_windows()
		ps1Filename = "$sourcePath/file/select_xml.ps1"
		str = readstring(`Powershell.exe -executionpolicy remotesigned -File $ps1Filename`)
	else
		print(message)
		str = readline()
	end
	xmlFilename = abspath(String(strip(chomp(str), [' ', '"']))) #"
	
	# check selected file
	@assert(ispath(xmlFilename), "file does not exist: '$xmlFilename'")
	@assert(split(xmlFilename, '.')[end] == "xml", "file is not xml file: '$xmlFilename'")
	
	return xmlFilename
end
