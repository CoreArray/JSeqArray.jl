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
	seqFilterPush, seqFilterPop, seqFilterGet, seqGetData, seqApply



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



####  Internal functions  ####

function gds_dim(file::TypeSeqArray)
	return ccall((:SEQ_GetSpace, LibSeqArray), Vector{Int64}, (Cint,),
		file.gds.id)
end

function gds_seldim(file::TypeSeqArray)
	return ccall((:SEQ_GetSelSpace, LibSeqArray), Vector{Int64}, (Cint,),
		file.gds.id)
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
function seqFilterSet(file::TypeSeqArray,
		sample_id::Union{Void, Vector} = nothing,
		variant_id::Union{Void, Vector} = nothing,
		intersect::Bool=false, verbose::Bool=true)
	# set samples
	if sample_id != nothing
		sampset = Set(sample_id)
		if !intersect
			seqFilterReset(file, true, false, false)
		end
		sampid = seqGetData(file, "sample.id")
		sampsel = [ in(x, sampset) for x in sampid ]
		seqFilterSet2(file, sampsel, nothing, intersect, verbose)
	end
	# set variants
	if variant_id != nothing
		varset = Set(variant_id)
		if !intersect
			seqFilterReset(file, false, true, false)
		end
		varid = seqGetData(file, "variant.id")
		varsel = [ in(x, varset) for x in varid ]
		seqFilterSet2(file, nothing, varsel, intersect, verbose)
	end
	return nothing
end



# Set a filter on variants or samples using an index vector or a logical vector
function seqFilterSet2(file::TypeSeqArray,
		sample::Union{Void, Vector{Bool}, Vector{Int}, UnitRange{Int}}=nothing,
		variant::Union{Void, Vector{Bool}, Vector{Int}, UnitRange{Int}}=nothing,
		intersect::Bool=false, verbose::Bool=true)
	# set samples
	if sample != nothing
		if typeof(sample) != Vector{Bool}
			if intersect
				flag = zeros(Bool, gds_seldim(file)[2])
			else
				flag = zeros(Bool, gds_dim(file)[2])
			end
			flag[sample] = true
			sample = flag
		end
		ccall((:SEQ_SetSampleB, LibSeqArray), Void,
			(Cint,Any,Bool,Bool), file.gds.id, sample, intersect, verbose)
	end
	# set variants
	if variant != nothing
		if typeof(variant) != Vector{Bool}
			if intersect
				flag = zeros(Bool, gds_seldim(file)[3])
			else
				flag = zeros(Bool, gds_dim(file)[3])
			end
			flag[variant] = true
			variant = flag
		end
		ccall((:SEQ_SetVariantB, LibSeqArray), Void,
			(Cint,Any,Bool,Bool), file.gds.id, variant, intersect, verbose)
	end
	return nothing
end



# Reset the filter
function seqFilterReset(file::TypeSeqArray, sample::Bool=true,
		variant::Bool=true, verbose::Bool=true)
	# set samples
	if sample
		ccall((:SEQ_SetSampleB, LibSeqArray), Void,
			(Cint,Any,Bool,Bool), file.gds.id, nothing, false, verbose)
	end
	# set variants
	if variant
		ccall((:SEQ_SetVariantB, LibSeqArray), Void,
			(Cint,Any,Bool,Bool), file.gds.id, nothing, false, verbose)
	end
	return nothing
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
	return ccall((:SEQ_GetFilter, LibSeqArray), Vector{Bool}, (Cint,Bool),
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
function seqApply(fun::Function, file::TypeSeqArray,
		name::Union{String, Vector{String}}, asis::String="none",
		verbose::Bool=false, bsize::Int=1024; args...)
	# build additional parameters for the user-defined function
	args = Vector{Any}([ x[2] for x in args ])
	# TODO: check the number of arguments
	# c call
	ans = ccall((:SEQ_BApply_Variant, LibSeqArray), Any,
		(Cint, Any, Function, Cstring, Cint, Bool, Vector{Any}),
		file.gds.id, name, fun, asis, bsize, verbose, args)
	if asis == "unlist"
		ans = vcat(ans...)
	end
	return ans
end


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
