
<a id='JSeqArray.jl-Documentation-1'></a>

# JSeqArray.jl Documentation

- [JSeqArray.jl Documentation](index.md#JSeqArray.jl-Documentation-1)
    - [Functions](index.md#Functions-1)




<a id='Functions-1'></a>

## Functions

<a id='JSeqArray.seqOpen' href='#JSeqArray.seqOpen'>#</a>
**`JSeqArray.seqOpen`** &mdash; *Function*.



```
seqOpen(filename, readonly, allow_dup)
```

Opens a SeqArray GDS file.

**Arguments**

  * `filename::String`: the file name of a SeqArray file
  * `readonly::Bool=true`: if true, the file is opened read-only; otherwise, it is allowed to write data to the file
  * `allow_dup::Bool=false`: if true, it is allowed to open a GDS file with read-only mode when it has been opened in the same session

**Examples**

```julia
julia> f = seqOpen(seqExample(:kg))

julia> f

julia> seqClose(f)
```

<a id='JSeqArray.seqClose-Tuple{JSeqArray.TypeSeqFile}' href='#JSeqArray.seqClose-Tuple{JSeqArray.TypeSeqFile}'>#</a>
**`JSeqArray.seqClose`** &mdash; *Method*.



```
seqClose(file)
```

Closes a SeqArray GDS file which is open.

**Arguments**

  * `file::TypeSeqFile`: a SeqArray julia object

<a id='JSeqArray.seqFilterSet-Tuple{JSeqArray.TypeSeqFile}' href='#JSeqArray.seqFilterSet-Tuple{JSeqArray.TypeSeqFile}'>#</a>
**`JSeqArray.seqFilterSet`** &mdash; *Method*.



```
seqFilterSet(file; sample_id, variant_id, intersect, verbose)
```

Sets a filter to sample and/or variant.

**Arguments**

  * `file::TypeSeqFile`: a SeqArray julia object
  * `sample_id::Union{Void, Vector}=nothing`: sample ID to be selected, or `nothing` for no action
  * `variant_id::Union{Void, Vector}=nothing`: variant ID to be selected, or `nothing` for no action
  * `intersect::Bool=false`: if false, the candidate samples/variants for selection are all samples/variants; if true, the candidate samples/variants are from the selected samples/variants defined via the previous call
  * `verbose::Bool=true`: if true, show information

**Examples**

```julia-repl
julia> f = seqOpen(seqExample(:kg));

julia> sid = seqGetData(f, "sample.id");

julia> vid = seqGetData(f, "variant.id");

julia> seqFilterSet(f, sample_id=sid[4:10], variant_id=vid[2:6])
Number of selected samples: 7
Number of selected variants: 5

julia> seqClose(f)
```

<a id='JSeqArray.seqFilterSet2-Tuple{JSeqArray.TypeSeqFile}' href='#JSeqArray.seqFilterSet2-Tuple{JSeqArray.TypeSeqFile}'>#</a>
**`JSeqArray.seqFilterSet2`** &mdash; *Method*.



```
seqFilterSet2(file; sample, variant, intersect, verbose)
```

Sets a filter to sample and/or variant.

**Arguments**

  * `file::TypeSeqFile`: a SeqArray julia object
  * `sample::Union{Void, Vector{Bool}, Vector{Int}, UnitRange{Int}}=nothing`: sample(s) to be selected, or `nothing` for no action
  * `variant::Union{Void, Vector{Bool}, Vector{Int}, UnitRange{Int}}=nothing`: variant(s) to be selected, or `nothing` for no action
  * `intersect::Bool=false`: if false, the candidate samples/variants for selection are all samples/variants; if true, the candidate samples/variants are from the selected samples/variants defined via the previous call
  * `verbose::Bool=true`: if true, show information

<a id='JSeqArray.seqFilterSplit-Tuple{JSeqArray.TypeSeqFile,Int64,Int64}' href='#JSeqArray.seqFilterSplit-Tuple{JSeqArray.TypeSeqFile,Int64,Int64}'>#</a>
**`JSeqArray.seqFilterSplit`** &mdash; *Method*.



```
seqFilterSplit(file, index, count; verbose)
```

Splits the variants into multiple parts equally and selects the specified part.

**Arguments**

  * `file::TypeSeqFile`: a SeqArray julia object
  * `index::Int`: selects the `index`th part (starting from 1)
  * `index::Int`: the total number of non-overlapping parts
  * `verbose::Bool=true`: if true, show information

**Details**

Users can define a subset of variants before calling `seqFilterSplit()` and split the selection of variants into multiple parts.

**Examples**

```julia-repl
julia> f = seqOpen(seqExample(:kg));

julia> seqFilterSplit(f, 2, 5)
Number of selected variants: 3,954

julia> seqClose(f)
```

<a id='JSeqArray.seqFilterReset-Tuple{JSeqArray.TypeSeqFile}' href='#JSeqArray.seqFilterReset-Tuple{JSeqArray.TypeSeqFile}'>#</a>
**`JSeqArray.seqFilterReset`** &mdash; *Method*.



```
seqFilterReset(file; sample, variant, verbose)
```

Resets the sample and variant filters.

**Arguments**

  * `file::TypeSeqFile`: a SeqArray julia object
  * `sample::Bool=true`: if true, resets the sample filter
  * `variant::Bool=true`: if true, resets the variant filter
  * `verbose::Bool=true`: if true, show information

<a id='JSeqArray.seqFilterPush' href='#JSeqArray.seqFilterPush'>#</a>
**`JSeqArray.seqFilterPush`** &mdash; *Function*.



```
seqFilterPush(file, reset)
```

Pushes the sample and variant filters to the stack for future uses.

**Arguments**

  * `file::TypeSeqFile`: a SeqArray julia object
  * `reset::Bool=false`: if true, reset the sample and variant filters

<a id='JSeqArray.seqFilterPop-Tuple{JSeqArray.TypeSeqFile}' href='#JSeqArray.seqFilterPop-Tuple{JSeqArray.TypeSeqFile}'>#</a>
**`JSeqArray.seqFilterPop`** &mdash; *Method*.



```
seqFilterPop(file)
```

Uses the last sample and variant filters saved in the stack, and removes them from the stack.

**Arguments**

  * `file::TypeSeqFile`: a SeqArray julia object

<a id='JSeqArray.seqFilterGet' href='#JSeqArray.seqFilterGet'>#</a>
**`JSeqArray.seqFilterGet`** &mdash; *Function*.



```
seqFilterGet(file, sample)
```

Gets the filter of samples and variants.

**Arguments**

  * `file::TypeSeqFile`: a SeqArray julia object
  * `sample::Bool=true`: if true, returns a logical vector for the sample filter; otherwise, returns a logical vector for the variant filter

<a id='JSeqArray.seqGetData-Tuple{JSeqArray.TypeSeqFile,String}' href='#JSeqArray.seqGetData-Tuple{JSeqArray.TypeSeqFile,String}'>#</a>
**`JSeqArray.seqGetData`** &mdash; *Method*.



```
seqGetData(file, name)
```

Gets data from a SeqArray GDS file.

**Arguments**

  * `file::TypeSeqFile`: a SeqArray julia object
  * `name::String`: the variable name, see the details

**Details**

The variable name should be

  * "sample.id", "variant.id", "position", "chromosome", "allele"
  * "genotype" for 3-dim UInt8 array (ploidy, sample, variant) where 0 is the reference allele, 1 is the first alternative allele, 0xFF is missing value
  * "annotation/id", "annotation/qual", "annotation/filter", "annotation/info/VARIABLE_NAME", "annotation/format/VARIABLE_NAME"
  * "#dosage" for a dosage matrix (sample, variant) of reference allele (UInt8: 0, 1 and 2 for diploid genotypes, 0xFF for missing values)
  * "#num_allele" returns an integer vector with the numbers of distinct alleles

**Examples**

```julia-repl
julia> f = seqOpen(seqExample(:kg));

julia> pos = seqGetData(f, "position"); println(typeof(pos), ", ", length(pos))
Array{Int32,1}, 19773

julia> geno = seqGetData(f, "genotype"); println(typeof(geno), ", ", size(geno))
Array{UInt8,3}, (2,1092,19773)

julia> dosage = seqGetData(f, "#dosage"); println(typeof(dosage), ", ", size(dosage))
Array{UInt8,2}, (1092,19773)

julia> seqClose(f)
```

<a id='JSeqArray.seqApply-Tuple{Function,JSeqArray.TypeSeqFile,Union{Array{String,1},String},Vararg{Any,N}}' href='#JSeqArray.seqApply-Tuple{Function,JSeqArray.TypeSeqFile,Union{Array{String,1},String},Vararg{Any,N}}'>#</a>
**`JSeqArray.seqApply`** &mdash; *Method*.



```
seqApply(fun, file, name, args...; asis, bsize, verbose, kwargs...)
```

Applies the user-defined function over array margins.

**Arguments**

  * `fun::Function`: the user-defined function
  * `file::TypeSeqFile`: a SeqArray julia object
  * `name::Union{String, Vector{String}}`: the variable name(s), see the details
  * `args`: the optional arguments passed to the user-defined function
  * `asis::Symbol=:none`: `:none` (no return), `:unlist` (returns a vector which contains all the atomic components) or `:list` (returns a vector according to each block)
  * `bsize::Int=1024`: block size for the number of variants in a block
  * `verbose::Bool=true`: if true, show progress information
  * `kwargs`: the keyword optional arguments passed to the user-defined function

**Details**

The variable name should be

  * "sample.id", "variant.id", "position", "chromosome", "allele"
  * "genotype" for 3-dim UInt8 array (ploidy, sample, variant) where 0 is the reference allele, 1 is the first alternative allele, 0xFF is missing value
  * "annotation/id", "annotation/qual", "annotation/filter", "annotation/info/VARIABLE_NAME", "annotation/format/VARIABLE_NAME"
  * "#dosage" for a dosage matrix (sample, variant) of reference allele (UInt8: 0, 1 and 2 for diploid genotypes, 0xFF for missing values)
  * "#num_allele" returns an integer vector with the numbers of distinct alleles

The algorithm is highly optimized by blocking the computations to exploit the high-speed memory instead of disk.

**Examples**

```julia
julia> f = seqOpen(seqExample(:kg))
julia> seqApply(f, "genotype", asis=:unlist) do geno
           return sum(geno)
       end
julia> seqClose(f)
```

<a id='JSeqArray.seqParallel-Tuple{Function,JSeqArray.TypeSeqFile,Vararg{Any,N}}' href='#JSeqArray.seqParallel-Tuple{Function,JSeqArray.TypeSeqFile,Vararg{Any,N}}'>#</a>
**`JSeqArray.seqParallel`** &mdash; *Method*.



```
seqParallel(fun, file, args...; split, combine, kwargs...)
```

Applies a user-defined function in parallel.

**Arguments**

  * `fun::Function`: the user-defined function
  * `file::TypeSeqFile`: a SeqArray julia object
  * `args`: the optional arguments passed to the user-defined function
  * `split::Symbol=:byvariant`: `:none` for no split, `:byvariant` for spliting the dataset by variant according to multiple processes
  * `combine::Union{Symbol, Function}=:unlist`: `:none` (no return), `:unlist` (returns a vector which contains all the atomic components) or `:list` (returns a vector according to each process)
  * `kwargs`: the keyword optional arguments passed to the user-defined function

**Details**

**Examples**

<a id='JSeqArray.seqAttr-Tuple{JSeqArray.TypeSeqFile,Symbol}' href='#JSeqArray.seqAttr-Tuple{JSeqArray.TypeSeqFile,Symbol}'>#</a>
**`JSeqArray.seqAttr`** &mdash; *Method*.



```
seqAttr(file, name)
```

Applies a user-defined function in parallel.

**Arguments**

  * `file::TypeSeqFile`: a SeqArray julia object
  * `name::Symbol`: the symbol name for a specified attribute

**Details**

`name::Symbol =`

  * `:nsampe` - the total number of samples
  * `:nselsamp` - the number of selected samples
  * `:nvar` - the total number of variants
  * `:nselvar` - the number of selected variants
  * `:ploidy` - the number of sets of chromosomes

**Examples**

```julia-repl
julia> f = seqOpen(seqExample(:kg));

julia> seqFilterSet2(f, sample=5:10, variant=31:40)
Number of selected samples: 6
Number of selected variants: 10

julia> seqAttr(f, :nsamp)
1092

julia> seqAttr(f, :nselsamp)
6

julia> seqAttr(f, :nvar)
19773

julia> seqAttr(f, :nselvar)
10

julia> seqAttr(f, :ploidy)
2

julia> seqClose(f)
```

<a id='JSeqArray.seqExample-Tuple{Symbol}' href='#JSeqArray.seqExample-Tuple{Symbol}'>#</a>
**`JSeqArray.seqExample`** &mdash; *Method*.



```
seqExample(file)
```

Returns the example SeqArray file.

**Arguments**

  * `file::Symbol`: specify which SeqArray file, it should be :kg for 1KG_phase1_release_v3_chr22.gds

**Examples**

```julia-repl
julia> fn = seqExample(:kg);

julia> basename(fn)
"1KG_phase1_release_v3_chr22.gds"
```

<a id='Base.show-Tuple{IO,JSeqArray.TypeSeqFile}' href='#Base.show-Tuple{IO,JSeqArray.TypeSeqFile}'>#</a>
**`Base.show`** &mdash; *Method*.



```
show(io, file; attr, all)
```

Applies a user-defined function in parallel.

**Arguments**

  * `io::`: I/O stream
  * `file::TypeSeqFile`: a SeqArray julia object
  * `attr::Bool=false`: if true, shows all attributes
  * `all::Bool=false`: if true, show all GDS nodes including hidden nodes

