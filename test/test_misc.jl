
## Test: the calculation of allele frequencies

f = seqOpen(seqExample(:kg))
println("Calculation of allele frequencies")

try
	af = seqApply(f, "genotype", asis=:unlist) do geno::Array{UInt8,3}
	    N = size(geno, 3)
    	rv = Vector{Float64}(N)
	    for k in 1:N
    	    sum = 0; n = 0
	        for g in geno[:,:,k]
    	        if g != 0xFF
        	        sum += g == 0  # 0 is the reference allele
            	    n += 1
	            end
	        end
	        rv[k] = sum / n
	    end
    	return rv
	end

	af_min = minimum(af)
	af_max = maximum(af)
	af_mean = mean(af)
	af_median = median(af)

	@test af_min ≈ 0.0
	@test af_max ≈ 1.0
	@test af_mean ≈ 0.9286052603644784
	@test af_median ≈ 0.9963369963369964

finally
	seqClose(f)
end



