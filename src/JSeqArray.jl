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

export TypeSeqFile, TypeVarData,
	seqExample, seqOpen, seqClose, seqFilterSet, seqFilterSet2, seqFilterSplit,
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
type TypeSeqFile <: anygdsfile
	gds::type_gdsfile
	auxiliary::Any
end

# Type for variable-length data
immutable TypeVarData
	index::Vector{Int32}
	data::Any
end



####  Internal functions  ####

# return (ploidy, total # of samples, total # of variants)
function gds_dim(file::TypeSeqFile)
	return ccall((:SEQ_GetSpace, LibSeqArray), Vector{Int64}, (Cint,),
		file.gds.id)
end

# return (ploidy, # of selected samples, # of selected variants)
function gds_seldim(file::TypeSeqFile)
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




####  Example GDS Files  ####

# Get the directory of jugds *.h header files
"""
	seqExample(file)
Returns the example SeqArray file.
# Arguments
* `file::Symbol`: specify which SeqArray file, it should be :kg for 1KG_phase1_release_v3_chr22.gds
# Examples
```jldoctest
julia> fn = seqExample(:kg);

julia> basename(fn)
"1KG_phase1_release_v3_chr22.gds"
```
"""
function seqExample(file::Symbol)
	path = joinpath(Pkg.dir(), "JSeqArray", "demo", "data")
	if file == :kg
		fn = joinpath(path, "1KG_phase1_release_v3_chr22.gds")
	else
		throw(ArgumentError("'file' should be :kg."))
	end
	return fn
end



####  GDS File  ####

# Open a SeqArray file
"""
	seqOpen(filename, readonly, allow_dup)
Opens a SeqArray GDS file.
# Arguments
* `filename::String`: the file name of a SeqArray file
* `readonly::Bool=true`: if true, the file is opened read-only; otherwise, it is allowed to write data to the file
* `allow_dup::Bool=false`: if true, it is allowed to open a GDS file with read-only mode when it has been opened in the same session
# Examples
```julia
julia> f = seqOpen(seqExample(:kg))

julia> f

julia> seqClose(f)
```
"""
function seqOpen(filename::String, readonly::Bool=true, allow_dup::Bool=false)
	ff = open_gds(filename, readonly, allow_dup)
	# TODO: check file structure
	ccall((:SEQ_File_Init, LibSeqArray), Void, (Cint,), ff.id)
	return TypeSeqFile(ff, nothing)
end



# Close the SeqArray file
"""
	seqClose(file)
Closes a SeqArray GDS file which is open.
# Arguments
* `file::TypeSeqFile`: a SeqArray julia object
"""
function seqClose(file::TypeSeqFile)
	fid = file.gds.id
	close_gds(file.gds)
	ccall((:SEQ_File_Done, LibSeqArray), Void, (Cint,), fid)
	return nothing
end



# Set a filter on variants or samples with sample or variant IDs
"""
	seqFilterSet(file; sample_id, variant_id, intersect, verbose)
Sets a filter to sample and/or variant.
# Arguments
* `file::TypeSeqFile`: a SeqArray julia object
* `sample_id::Union{Void, Vector}=nothing`: sample ID to be selected, or nothing for no action
* `variant_id::Union{Void, Vector}=nothing`: variant ID to be selected, or nothing for no action
* `intersect::Bool=false`: if false, the candidate samples/variants for selection are all samples/variants; if true, the candidate samples/variants are from the selected samples/variants defined via the previous call
* `verbose::Bool=true`: if true, show information
# Examples
```julia
julia> f = seqOpen(seqExample(:kg))
julia> sid = seqGetData(f, "sample.id")
julia> vid = seqGetData(f, "variant.id")
julia> seqFilterSet(f, sample_id=sid[4:10], variant_id=vid[2:6])
julia> seqClose(f)
```
"""
function seqFilterSet(file::TypeSeqFile;
		sample_id::Union{Void, Vector} = nothing,
		variant_id::Union{Void, Vector} = nothing,
		intersect::Bool=false, verbose::Bool=true)
	# set samples
	if sample_id != nothing
		sampset = Set(sample_id)
		if !intersect
			seqFilterReset(file, sample=true, variant=false, verbose=false)
		end
		sampid = seqGetData(file, "sample.id")
		sampsel = [ in(x, sampset) for x in sampid ]
		seqFilterSet2(file, sample=sampsel, intersect=intersect, verbose=verbose)
	end
	# set variants
	if variant_id != nothing
		varset = Set(variant_id)
		if !intersect
			seqFilterReset(file, sample=false, variant=true, verbose=false)
		end
		varid = seqGetData(file, "variant.id")
		varsel = [ in(x, varset) for x in varid ]
		seqFilterSet2(file, variant=varsel, intersect=intersect, verbose=verbose)
	end
	return nothing
end



# Set a filter on variants or samples using an index vector or a logical vector
"""
	seqFilterSet2(file; sample, variant, intersect, verbose)
Sets a filter to sample and/or variant.
# Arguments
* `file::TypeSeqFile`: a SeqArray julia object
* `sample::Union{Void, Vector{Bool}, Vector{Int}, UnitRange{Int}}=nothing`: sample(s) to be selected, or nothing for no action
* `variant::Union{Void, Vector{Bool}, Vector{Int}, UnitRange{Int}}=nothing`: variant(s) to be selected, or nothing for no action
* `intersect::Bool=false`: if false, the candidate samples/variants for selection are all samples/variants; if true, the candidate samples/variants are from the selected samples/variants defined via the previous call
* `verbose::Bool=true`: if true, show information
"""
function seqFilterSet2(file::TypeSeqFile;
		sample::Union{Void, Vector{Bool}, Vector{Int}, UnitRange{Int}}=nothing,
		variant::Union{Void, Vector{Bool}, Vector{Int}, UnitRange{Int}}=nothing,
		intersect::Bool=false, verbose::Bool=true)
	# set samples
	if sample != nothing
		if !isa(sample, Vector{Bool})
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
		if !isa(variant, Vector{Bool})
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
"""
	seqFilterSplit(file, index, count; verbose)
Splits the variants into multiple parts equally and selects the specified part.
# Arguments
* `file::TypeSeqFile`: a SeqArray julia object
* `index::Int`: selects the `index`th part (starting from 1)
* `index::Int`: the total number of non-overlapping parts
* `verbose::Bool=true`: if true, show information
# Details
Users can define a subset of variants before calling `seqFilterSplit()` and split the selection of variants into multiple parts.
"""
function seqFilterSplit(file::TypeSeqFile, index::Int, count::Int;
		verbose::Bool=true)
	if count < 1
		throw(ArgumentError("'count' should be > 0."))
	end
	if index < 1 || index > count
		throw(ArgumentError("'index' should be between 1 and $count."))
	end
	ss = split_count(gds_seldim(file)[3], count)
	seqFilterSet2(file, variant=ss[index], intersect=true, verbose=verbose)
	return nothing
end



# Reset the filter
"""
	seqFilterReset(file; sample, variant, verbose)
Resets the sample and variant filters.
# Arguments
* `file::TypeSeqFile`: a SeqArray julia object
* `sample::Bool=true`: if true, resets the sample filter
* `variant::Bool=true`: if true, resets the variant filter
* `verbose::Bool=true`: if true, show information
"""
function seqFilterReset(file::TypeSeqFile; sample::Bool=true,
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
"""
	seqFilterPush(file, reset)
Pushes the sample and variant filters to the stack for future uses.
# Arguments
* `file::TypeSeqFile`: a SeqArray julia object
* `reset::Bool=false`: if true, reset the sample and variant filters
"""
function seqFilterPush(file::TypeSeqFile, reset::Bool=false)
	ccall((:SEQ_FilterPush, LibSeqArray), Void, (Cint,Bool), file.gds.id, reset)
	return nothing
end



# Pop a filter
"""
	seqFilterPop(file)
Uses the last sample and variant filters saved in the stack, and removes them from the stack.
# Arguments
* `file::TypeSeqFile`: a SeqArray julia object
"""
function seqFilterPop(file::TypeSeqFile)
	ccall((:SEQ_FilterPop, LibSeqArray), Void, (Cint,), file.gds.id)
	return nothing
end



# Get a sample/variant filter
"""
	seqFilterGet(file, sample)
Gets the filter of samples and variants.
# Arguments
* `file::TypeSeqFile`: a SeqArray julia object
* `sample::Bool=true`: if true, returns a logical vector for the sample filter; otherwise, returns a logical vector for the variant filter
"""
function seqFilterGet(file::TypeSeqFile, sample::Bool=true)
	return ccall((:SEQ_GetFilter, LibSeqArray), Vector{Bool}, (Cint,Bool),
		file.gds.id, sample)
end



# Get data
"""
	seqGetData(file, name)
Gets data from a SeqArray GDS file.
# Arguments
* `file::TypeSeqFile`: a SeqArray julia object
* `name::String`: the variable name, see the details
# Details
The variable name should be
* "sample.id", "variant.id", "position", "chromosome", "allele"
* "genotype" for 3-dim UInt8 array (ploidy, sample, variant) where 0 is the reference allele, 1 is the first alternative allele, 0xFF is missing value
* "annotation/id", "annotation/qual", "annotation/filter", "annotation/info/VARIABLE_NAME", "annotation/format/VARIABLE_NAME"
* "#dosage" for a dosage matrix (sample, variant) of reference allele (UInt8: 0, 1 and 2 for diploid genotypes, 0xFF for missing values)
* "#num_allele" returns an integer vector with the numbers of distinct alleles
# Examples
```jldoctest
julia> f = seqOpen(seqExample(:kg));

julia> pos = seqGetData(f, "position"); println(typeof(pos), ", ", length(pos))
Array{Int32,1}, 19773

julia> geno = seqGetData(f, "genotype"); println(typeof(geno), ", ", size(geno))
Array{UInt8,3}, (2,1092,19773)

julia> dosage = seqGetData(f, "#dosage"); println(typeof(dosage), ", ", size(dosage))
Array{UInt8,2}, (1092,19773)

julia> seqClose(f)
```
"""
function seqGetData(file::TypeSeqFile, name::String)
	rv = ccall((:SEQ_GetData, LibSeqArray), Any, (Cint,Cstring),
		file.gds.id, name)
	if isa(rv, Vector{Any})
		rv = TypeVarData(rv[1], rv[2])
	end
	return rv
end



# Apply function over array margins
"""
	seqApply(fun, file, name, args...; asis, bsize, verbose, kwargs...)
Applies the user-defined function over array margins.
# Arguments
* `fun::Function`: the user-defined function
* `file::TypeSeqFile`: a SeqArray julia object
* `name::Union{String, Vector{String}}`: the variable name(s), see the details
* `args`: the optional arguments passed to the user-defined function
* `asis::Symbol=:none`: `:none` (no return), `:unlist` (returns a vector which contains all the atomic components) or `:list` (returns a vector according to each block)
* `bsize::Int=1024`: block size for the number of variants in a block
* `verbose::Bool=true`: if true, show progress information
* `kwargs`: the keyword optional arguments passed to the user-defined function
# Details
The variable name should be
* "sample.id", "variant.id", "position", "chromosome", "allele"
* "genotype" for 3-dim UInt8 array (ploidy, sample, variant) where 0 is the reference allele, 1 is the first alternative allele, 0xFF is missing value
* "annotation/id", "annotation/qual", "annotation/filter", "annotation/info/VARIABLE_NAME", "annotation/format/VARIABLE_NAME"
* "#dosage" for a dosage matrix (sample, variant) of reference allele (UInt8: 0, 1 and 2 for diploid genotypes, 0xFF for missing values)
* "#num_allele" returns an integer vector with the numbers of distinct alleles
The algorithm is highly optimized by blocking the computations to exploit the high-speed memory instead of disk.
# Examples
```julia
julia> f = seqOpen(seqExample(:kg))
julia> seqApply(f, "genotype", asis=:unlist) do geno
           return sum(geno)
       end
julia> seqClose(f)
```
"""
function seqApply(fun::Function, file::TypeSeqFile,
		name::Union{String, Vector{String}}, args...; asis::Symbol=:none,
		bsize::Int=1024, verbose::Bool=true, kwargs...)
	# check
	if asis!=:none && asis!=:unlist && asis!=:list
		throw(ArgumentError("'asis' should be :none, :unlist or :list."))
	end
	if bsize <= 0
		throw(ArgumentError("'bsize' should be greater than 0."))
	end
	if isa(name, String)
		name = [ name ]
	end
	# initialize
	dm = gds_seldim(file)[3]  # the number of selected variants
	bnum = div(dm, bsize)
	bnum += mod(dm, bsize) != 0
	rv = asis==:none ? nothing : Vector{Any}(bnum)
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
			seqFilterSet2(file, variant=idx[st:ed], verbose=false)
			st += bsize
			x = [ seqGetData(file, nm) for nm in name ]
			v = fun(x..., args...; kwargs...)
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
	if asis == :unlist
		rv = vcat(rv...)
	end
	return rv
end



####  Parallel functions  ####

# internal variables used for identifying processes
process_index = 0
process_count = 0


# Apply Functions in Parallel
"""
	seqParallel(fun, file, args...; split, combine, kwargs...)
Applies a user-defined function in parallel.
# Arguments
* `fun::Function`: the user-defined function
* `file::TypeSeqFile`: a SeqArray julia object
* `args`: the optional arguments passed to the user-defined function
* `split::Symbol=:byvariant`: `:none` for no split, `:byvariant` for spliting the dataset by variant according to multiple processes
* `combine::Union{Symbol, Function}=:unlist`: `:none` (no return), `:unlist` (returns a vector which contains all the atomic components) or `:list` (returns a vector according to each process)
* `kwargs`: the keyword optional arguments passed to the user-defined function
# Details
# Examples
"""
function seqParallel(fun::Function, file::TypeSeqFile, args...;
		split::Symbol=:byvariant, combine::Union{Symbol, Function}=:unlist,
		kwargs...)
	# check
	if split!=:byvariant && split!=:none
		throw(ArgumentError("'split' should be :byvaraint or :none."))
	end
	if isa(combine, Symbol)
		if combine!=:none && combine!=:unlist && combine!=:list
			throw(ArgumentError("'combine' should be :none, :unlist or :list."))
		end
	end
	# set remotecall
	@everywhere using JSeqArray
	ws = workers()
	rc = Vector{Any}(length(ws))
	fn = file.gds.filename
	ssel = seqFilterGet(file, true)
	vsel = seqFilterGet(file, false)
	for i in 1:length(ws)
		rc[i] = remotecall(ws[i], i, length(ws), fn, ssel, vsel, fun, split,
					args, kwargs) do i, cnt, fn, ssel, vsel, fun, split, args, kwargs
			process_index = i
			process_count = cnt
			rv = nothing
			gdsfile = seqOpen(fn, true, true)
			seqFilterSet2(gdsfile, sample=ssel, variant=vsel, verbose=false)
			if split==:byvariant
				seqFilterSplit(gdsfile, i, cnt, verbose=false)
			end
			try
				rv = fun(gdsfile, args...; kwargs...)
			finally
				seqClose(gdsfile)
			end
			return rv
		end
	end
	# remote run
	if isa(combine, Symbol)
		rv = [ fetch(r) for r in rc ]
		err = false
		for i in rv
			if isa(i, RemoteException)
				show(i)
				err = true
			end
		end
		if err
			error("RemoteException")
		end
		rv = combine==:none ? nothing : (combine==:unlist ? vcat(rv...) : rv)
	else
		rv = reduce(combine, map(fetch, rc))
	end
	# output
	return rv
end




####  Display  ####

function show(io::IO, file::TypeSeqFile; attr=false, all=false)
	print_with_color(:bold, io, "SeqArray ")
	show(io, file.gds, attr=attr, all=all)
end


end
