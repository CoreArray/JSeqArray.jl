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
	seqOpen, seqClose, seqFilterSet, seqFilterSet2, seqFilterSplit,
	seqFilterReset, seqFilterPush, seqFilterPop, seqFilterGet, seqGetData,
	seqApply, seqParallel



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

# ploidy X total sample X total variant
function gds_dim(file::TypeSeqArray)
	return ccall((:SEQ_GetSpace, LibSeqArray), Vector{Int64}, (Cint,),
		file.gds.id)
end

# ploidy X selected sample X selected variant
function gds_seldim(file::TypeSeqArray)
	return ccall((:SEQ_GetSelSpace, LibSeqArray), Vector{Int64}, (Cint,),
		file.gds.id)
end

# split total count
function split_count(total_count::Int64, count::Int)
	scale = total_count / count
	start = 1.0
	rv = Vector{UnitRange{Int}}(count)
	for i in 1:count
		st = round(start)
		start += scale
		ed = round(start) - 1
		rv[i] = UnitRange{Int}(st, ed)
	end
	return rv
end



####  Internal Progress Bar  ####

# progress bar
type progress_bar
	obj::Ptr{Void}
end

function progress_init(count::Int64, verbose::Bool)
	ptr = ccall((:SEQ_ProgressInit, LibSeqArray), Ptr{Void}, (Int64, Bool),
		count, verbose)
	return progress_bar(ptr)
end

function progress_forward(bar::progress_bar)
	ccall((:SEQ_ProgressForward, LibSeqArray), Void, (Ptr{Void},), bar.obj)
	return nothing
end

function progress_done(bar::progress_bar)
	ccall((:SEQ_ProgressDone, LibSeqArray), Void, (Ptr{Void},), bar.obj)
	bar.obj = C_NULL
	return nothing
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
	fid = file.gds.id
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
		ccall((:SEQ_SetSample, LibSeqArray), Void,
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
		ccall((:SEQ_SetVariant, LibSeqArray), Void,
			(Cint,Any,Bool,Bool), file.gds.id, variant, intersect, verbose)
	end
	return nothing
end



# Set a filter on variants within a block given by the total number of blocks
function seqFilterSplit(file::TypeSeqArray, index::Int, count::Int,
		verbose::Bool=true)
	if count < 1
		error("'count' should be > 0.")
	end
	if index < 1 || index > count
		error("'index' should be between 1 and $count.")
	end
	ss = split_count(gds_seldim(file)[3], count)
	seqFilterSet2(file, nothing, ss[index], true, verbose)
	return nothing
end



# Reset the filter
function seqFilterReset(file::TypeSeqArray, sample::Bool=true,
		variant::Bool=true, verbose::Bool=true)
	# set samples
	if sample
		ccall((:SEQ_SetSample, LibSeqArray), Void,
			(Cint,Any,Bool,Bool), file.gds.id, nothing, false, verbose)
	end
	# set variants
	if variant
		ccall((:SEQ_SetVariant, LibSeqArray), Void,
			(Cint,Any,Bool,Bool), file.gds.id, nothing, false, verbose)
	end
	return nothing
end



# Push a filter
function seqFilterPush(file::TypeSeqArray, reset::Bool=false)
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
	rv = ccall((:SEQ_GetData, LibSeqArray), Any, (Cint,Cstring),
		file.gds.id, name)
	if typeof(rv) == Vector{Any}
		rv = TypeVarData(rv[1], rv[2])
	end
	return rv
end



# Apply function over array margins
function seqApply(fun::Function, file::TypeSeqArray,
		name::Union{String, Vector{String}}; asis::String="none",
		verbose::Bool=true, bsize::Int=1024, args...)
	# check
	if asis!="none" && asis!="unlist" && asis!="list"
		error("'asis' should be \"none\", \"unlist\" or \"list\".")
	end
	if bsize <= 0
		error("'bsize' should be greater than 0.")
	end
	if typeof(name) == String
		name = [ name ]
	end
	# initialize
	dm = gds_seldim(file)[3]  # the number of selected variants
	bnum = div(dm, bsize)
	bnum += mod(dm, bsize) != 0
	rv = asis=="none" ? nothing : Vector{Any}(bnum)
	# build additional parameters for the user-defined function
	args = Vector{Any}([ x[2] for x in args ])
	# run
	seqFilterPush(file)
	progress = progress_init(bnum, verbose)
	try
		idx = find(seqFilterGet(file, false))
		st = 1
		
		for i in 1:bnum
			ed = st + bsize - 1
			if ed > dm
				ed = dm
			end
			seqFilterSet2(file, nothing, idx[st:ed], false, false)
			st += bsize
			x = [ seqGetData(file, nm) for nm in name ]
			v = fun(x..., args...)
			if rv != nothing
				rv[i] = v
			end
			progress_forward(progress)
		end
	finally
		seqFilterPop(file)
		progress_done(progress)
	end
	# output
	if asis == "unlist"
		rv = vcat(rv...)
	end
	return rv
end



####  Parallel functions  ####

# internal variables used for identifying processes
process_index = 0
process_count = 0


# Apply Functions in Parallel
function seqParallel(fun::Function, file::TypeSeqArray;
		split::String="by.variant", asis::Union{String, Function}="unlist",
		args...)
	# check
	if split!="by.variant" && split!="none"
		error("'split' should be \"by.varaint\" or \"none\".");
	end
	if typeof(asis) == String
		if asis!="none" && asis!="unlist" && asis!="list"
			error("'asis' should be \"none\", \"unlist\" or \"list\".")
		end
	end
	# set remotecall
	@everywhere using JSeqArray
	ws = workers()
	rc = Vector{Any}(length(ws))
	fn = file.gds.filename
	sel = seqFilterGet(file, false)
	for i in 1:length(ws)
		rc[i] = remotecall(ws[i], i, length(ws), fn, sel, fun, split, args) do i, cnt, fn, sel, fun, split, args
			process_index = i
			process_count = cnt
			rv = nothing
			gdsfile = seqOpen(fn, true, true)
			seqFilterSet2(gdsfile, nothing, sel, false, false)
			if split=="by.variant"
				seqFilterSplit(gdsfile, i, cnt, false)
			end
			try
				rv = fun(gdsfile, args...)
			finally
				seqClose(gdsfile)
			end
			return rv
		end
	end
	# remote run
	if typeof(asis) == String
		rv = [ fetch(r) for r in rc ]
		err = false
		for i in rv
			if typeof(i) == RemoteException
				println(i)
				err = true
			end
		end
		if err
			error("RemoteException")
		end
		rv = asis=="none" ? nothing : (asis=="unlist" ? vcat(rv...) : rv)
	else
		rv = reduce(asis, map(fetch, rc))
	end
	# output
	return rv
end




####  Display  ####

function show(io::IO, file::TypeSeqArray; attr=false, all=false)
	print_with_color(:bold, io, "SeqArray ")
	show(io, file.gds, attr=attr, all=all)
end


end
