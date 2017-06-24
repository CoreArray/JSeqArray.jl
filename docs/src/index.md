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
seqClose(file::TSeqGDSFile)
```

```@docs
seqFilterSet(file::TSeqGDSFile; sample_id::Union{Void, Vector} = nothing, variant_id::Union{Void, Vector} = nothing, intersect::Bool=false, verbose::Bool=true)
```

```@docs
seqFilterSet2(file::TSeqGDSFile; sample::Union{Void, Vector{Bool}, Vector{Int}, UnitRange{Int}}=nothing, variant::Union{Void, Vector{Bool}, Vector{Int}, UnitRange{Int}}=nothing, intersect::Bool=false, verbose::Bool=true)
```

```@docs
seqFilterSplit(file::TSeqGDSFile, index::Int, count::Int; verbose::Bool=true)
```

```@docs
seqFilterReset(file::TSeqGDSFile; sample::Bool=true, variant::Bool=true, verbose::Bool=true)
```

```@docs
seqFilterPush(file::TSeqGDSFile, reset::Bool=false)
```

```@docs
seqFilterPop(file::TSeqGDSFile)
```

```@docs
seqFilterGet(file::TSeqGDSFile, sample::Bool=true)
```

```@docs
seqGetData(file::TSeqGDSFile, name::String)
```

```@docs
seqApply(fun::Function, file::TSeqGDSFile, name::Union{String, Vector{String}}, args...; asis::Symbol=:none, bsize::Int=1024, verbose::Bool=true, kwargs...)
```

```@docs
seqParallel(fun::Function, file::TSeqGDSFile, args...; split::Symbol=:byvariant, combine::Union{Symbol, Function}=:unlist, kwargs...)
```

```@docs
seqAttr(file::TSeqGDSFile, name::Symbol)
```

```@docs
seqExample(file::Symbol)
```

```@docs
show(io::IO, file::TSeqGDSFile; attr=false, all=false)
```
