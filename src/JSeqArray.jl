# ===========================================================================
#
# JSeqArray.jl: Julia Interface to SeqArray Files
#
# Copyright (C) 2017    Xiuwen Zheng
#
# This is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License Version 3 as
# published by the Free Software Foundation.
#
# JSeqArray is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public
# License along with JSeqArray.
# If not, see <http://www.gnu.org/licenses/>.


module JSeqArray

using jugds

import Base: joinpath, show, print_with_color, println
import jugds: type_gdsfile, open_gds, close_gds, show

export TypeSeqArray, TypeVarData,
	seqOpen, seqClose, seqFilterSet, seqFilterSet2, seqFilterReset,
	seqFilterPush, seqFilterPop, seqFilterGet, seqGetData



####  Open and initialize the SeqArray binary library  ####

@static if is_apple()
	libfn = "libJSeqArray.dylib"
elseif is_windows()
	libfn = "libJSeqArray.dll"
elseif is_unix()
	libfn = "libJSeqArray.so"
else
	error("The platform is not supported.")
end

global libname = joinpath(Pkg.dir(), "JSeqArray", "deps", libfn)

if !isfile(libname)
	error("The SeqArray library cannot be found, please try Pkg.build(\"JSeqArray\").")
end

const LibSeqArray = libname

function __init__()
	ccall((:Init_GDS_Routines, LibSeqArray), Void, (Ptr{Void},), jugds.lib_c_api)
end



####  Type of GDS File and Node	 ####

# Type for a SeqArray file
type TypeSeqArray <: anygdsfile
	gds::type_gdsfile
	auxiliary::Any
end

# Type for variable-length data
type TypeVarData
	index::Vector{Int32}
	data::Any
end



####  GDS File  ####

# Open a SeqArray file
function seqOpen(filename::String, readonly::Bool=true, allow_dup::Bool=false)
	ff = open_gds(filename, readonly, allow_dup)
	# TODO: check file structure
	ccall((:SEQ_File_Init, LibSeqArray), Void, (Cint,), ff.id)
	return TypeSeqArray(ff, nothing)
end


# Close the SeqArray file
function seqClose(file::TypeSeqArray)
	fid = file.id
	close_gds(file.gds)
	ccall((:SEQ_File_Done, LibSeqArray), Void, (Cint,), fid)
	return nothing
end


# Set a filter on variants or samples with sample or variant IDs
function seqFilterSet(file::TypeSeqArray, sample_id=nothing, variant_id=nothing,
	intersect::Bool=false, verbose::Bool=true)

end


# Set a filter on variants or samples with an index vector
function seqFilterSet2(file::TypeSeqArray, sample_id=nothing, variant_id=nothing,
	intersect::Bool=false, verbose::Bool=true)

end


# Reset the filter
function seqFilterReset(file::TypeSeqArray, sample::Bool=true,
	variant::Bool=true, verbose::Bool=true)

end


# Push a filter
function seqFilterPush(file::TypeSeqArray, reset::Bool=true)
	ccall((:SEQ_FilterPush, LibSeqArray), Void, (Cint,Bool), file.gds.id, reset)
	return nothing
end


# Pop a filter
function seqFilterPop(file::TypeSeqArray)
	ccall((:SEQ_FilterPop, LibSeqArray), Void, (Cint,), file.gds.id)
	return nothing
end


# Get a sample/variant filter
function seqFilterGet(file::TypeSeqArray, sample::Bool=true)
	return ccall((:SEQ_GetSpace, LibSeqArray), Vector{Bool}, (Cint,Bool),
		file.gds.id, sample)
end

# Get data
function seqGetData(file::TypeSeqArray, name::String)
	ans = ccall((:SEQ_GetData, LibSeqArray), Any, (Cint,Cstring),
		file.gds.id, name)
	if typeof(ans) == Vector{Any}
		ans = TypeVarData(ans[1], ans[2])
	end
	return ans
end


# Apply function over array margins
# function seqApply(file::TypeSeqArray, name::String, fun, param=nothing,
# 	as_is='none', bsize=1024, verbose::Bool=false)
# end


# Apply Functions in Parallel
# function seqRunParallel(file::TypeSeqArray, fun, param=nothing, ncpu=0,
#	split='by.variant', combine='unlist')
# end




####  Display  ####

function show(io::IO, file::TypeSeqArray, attr=false, all=false)
	print_with_color(:bold, io, "SeqArray ")
	show(io, file.gds, attr, all)
end


end
