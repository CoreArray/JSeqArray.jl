# JSeqArray.jl Documentation

```@contents
```

```@meta
DocTestSetup = quote
	using JSeqArray
end
```


## Functions

```@docs
seqOpen(filename::String, readonly::Bool=true, allow_dup::Bool=false)
```

```@docs
seqClose(file::TypeSeqFile)
```

```@docs
seqFilterSet(file::TypeSeqFile; sample_id::Union{Void, Vector} = nothing, variant_id::Union{Void, Vector} = nothing, intersect::Bool=false, verbose::Bool=true)
```

```@docs
seqFilterSet2(file::TypeSeqFile; sample::Union{Void, Vector{Bool}, Vector{Int}, UnitRange{Int}}=nothing, variant::Union{Void, Vector{Bool}, Vector{Int}, UnitRange{Int}}=nothing, intersect::Bool=false, verbose::Bool=true)
```

```@docs
seqFilterSplit(file::TypeSeqFile, index::Int, count::Int; verbose::Bool=true)
```

```@docs
seqFilterReset(file::TypeSeqFile; sample::Bool=true, variant::Bool=true, verbose::Bool=true)
```

```@docs
seqFilterPush(file::TypeSeqFile, reset::Bool=false)
```

```@docs
seqFilterPop(file::TypeSeqFile)
```

```@docs
seqFilterGet(file::TypeSeqFile, sample::Bool=true)
```

```@docs
seqGetData(file::TypeSeqFile, name::String)
```

```@docs
seqApply(fun::Function, file::TypeSeqFile, name::Union{String, Vector{String}}, args...; asis::Symbol=:none, bsize::Int=1024, verbose::Bool=true, kwargs...)
```

```@docs
seqParallel(fun::Function, file::TypeSeqFile, args...; split::Symbol=:byvariant, combine::Union{Symbol, Function}=:unlist, kwargs...)
```

```@docs
seqExample(file::Symbol)
```

```@docs
show(io::IO, file::TypeSeqFile; attr=false, all=false)
```
