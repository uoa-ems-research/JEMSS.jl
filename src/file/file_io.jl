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

# misc file input/output functions

Base.adjoint(s::T) where {T<:AbstractString} = s # not sure where else to put this

function readDlmFileNextLine!(file::IOStream; delim::Char=delimiter)
    try
        return readdlm(IOBuffer(readline(file)), delim)
    catch
        return []
    end
end

function readDlmFile(filename::String; delim::Char=delimiter)
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
function writeDlmLine!(file::IOStream, args...; delim::Char=delimiter)
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

struct Table
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
        @assert(length(header) == numCols)
        @assert(allunique(header))

        headerDict = arrayDict(header)
        columns = Dict{Any,Vector{Any}}()
        for i = 1:numCols
            columns[header[i]] = data[:, i]
            # columns[header[i]] = view(data, :, i)
        end

        return new(name, header, headerDict, data, columns)
    end

    function Table(name, header; rows=[], cols=[])
        @assert(isempty(rows) + isempty(cols) == 1)
        if !isempty(rows)
            return Table(name, header, adjoint(hcat(rows...)))
        elseif !isempty(cols)
            return Table(name, header, hcat(cols...))
        end
    end
end

function writeTablesToFile!(file::IOStream, table::Table;
    writeNumRows::Bool=false, writeNumCols::Bool=false)

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

function writeTablesToFile!(file::IOStream, tables::Dict{String,Table})
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

function writeTablesToFile(filename::String, tables::Dict{String,Table})
    file = open(filename, "w")
    writeTablesToFile!(file, tables)
    close(file)
end

function readTablesFromFile(filename::String)
    return readTablesFromData(readDlmFile(filename))
end

function readTablesFromData(data::Array{Any,2})
    (numDataRows, numDataCols) = size(data)

    tables = Dict{String,Table}()
    i = 1
    while i <= numDataRows
        # find start of table
        while i <= numDataRows && isempty(data[i, 1])
            i += 1
        end
        if i > numDataRows
            break
        end

        tableName = data[i, 1]

        # read row and column counts, if given
        numRows = numCols = nullIndex
        (size(data, 2) >= 2) && !isempty(data[i, 2]) && (numRows = data[i, 2])
        (size(data, 2) >= 3) && !isempty(data[i, 3]) && (numCols = data[i, 3])

        i += 1 # move to table header

        # determine column count
        j = 1
        while j <= numDataCols && !isempty(data[i, j])
            j += 1
        end
        numCols == nullIndex ? numCols = j - 1 : @assert(numCols == j - 1)

        # read table header
        tableHeader = data[i, 1:numCols]

        i += 1 # move to table data

        # determine row count
        j = i
        while j <= numDataRows && !isempty(data[j, 1])
            j += 1
        end
        numRows == nullIndex ? numRows = j - i : @assert(numRows == j - i)

        # create table
        tableData = data[(i-1).+(1:numRows), 1:numCols]
        table = Table(convert(String, tableName), tableHeader, tableData)

        # add table to tables
        @assert(!haskey(tables, table.name))
        tables[table.name] = table

        # go to next table
        i += numRows # i = j
        @assert(i > numDataRows || isempty(data[i, 1]))
    end

    return tables
end

# given a table and names of fields from table header to use,
# return a vector of dicts for each row, with each dict mapping from
# a field name to the value for that row and field
function tableRowsFieldDicts(table::Table, fieldNames::Vector{T}) where {T}
    fieldNames = convert(Vector{String}, fieldNames)
    fieldsDict = Dict{String,Any}([f => nothing for f in fieldNames])
    numRows = size(table.data, 1)
    rowsFieldDicts = [deepcopy(fieldsDict) for i = 1:numRows] # rowsFieldDicts[i] = dict with fieldName => table.data[i,table.headerDict[fieldName]] for fieldName in fieldNames
    for fieldName in fieldNames
        j = table.headerDict[fieldName]
        for i = 1:numRows
            rowsFieldDicts[i][fieldName] = table.data[i, j]
        end
    end
    return rowsFieldDicts
end

# parse json attributes column from table
function parseAttributesColumn(table::Table)
    if haskey(table.columns, "attributes")
        return [isempty(s) ? Dict{String,Any}() : JSON.parse(s) for s in table.columns["attributes"]]
    else
        return [Dict{String,Any}() for i = 1:size(table.data, 1)]
    end
end

# crc checksum
function fileChecksum(filename::String)
    return CRC32c.crc32c(Mmap.mmap(filename))
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

# if path is absolute then return it, otherwise prepend with (assumed) absolute path
joinPathIfNotAbs(absPath::String, path::String) = isabspath(path) ? abspath(path) : joinpath(absPath, path)

function interpolateString(s::String)
    return string("\"", escape_string(s), "\"") |> Meta.parse |> eval # note that eval only works for global vars
end
interpolateString(s::SubString{String}) = interpolateString(String(s))

# some convenient functions for reading xml files
xmlFileRoot(filename::String) = root(parse_file(interpolateString(filename)))
findElt = find_element # shorthand
xName(x) = LightXML.name(x) # to avoid conflict with JuMP.name
containsElt(parentElt::XMLElement, eltString::String) = findElt(parentElt, eltString) !== nothing
eltContent(elt::XMLElement) = content(elt)
eltContentVal(elt::XMLElement) = eval(Meta.parse(eltContent(elt)))
eltContentInterpVal(elt::XMLElement) = interpolateString(eltContent(elt))
eltContent(parentElt::XMLElement, eltString::String) =
    try
        content(findElt(parentElt, eltString))
    catch
        error("Element not found: $eltString")
    end
eltContentVal(parentElt::XMLElement, eltString::String) = eval(Meta.parse(eltContent(parentElt, eltString)))
eltContentInterpVal(parentElt::XMLElement, eltString::String) = interpolateString(eltContent(parentElt, eltString))
eltAttr = attribute
eltAttrVal(elt::XMLElement, attrString::String) = eval(Meta.parse(eltAttr(elt, attrString)))

function childrenNodeNames(parentElt::XMLElement)
    childNodes = Vector{String}()
    for childNode in child_nodes(parentElt)
        if is_elementnode(childNode)
            push!(childNodes, xName(childNode))
        end
    end
    return childNodes
end

function selectXmlFile(; message::String="Enter xml filename: ")
    if Sys.iswindows()
        ps1Filename = "$sourceDir/file/select_xml.ps1"
        str = read(`Powershell.exe -executionpolicy remotesigned -File $ps1Filename`, String)
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
