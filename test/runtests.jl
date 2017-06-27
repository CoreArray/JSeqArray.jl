using Base.Test
using JSeqArray
using StatsBase

tests = [ "misc" ]

for t in tests
    fn = joinpath(dirname(@__FILE__), "test_$t.jl")
    println("test/test_$t.jl ...")
    include(fn)
end


