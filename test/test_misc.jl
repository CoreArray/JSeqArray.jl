
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

	s = summarystats(af)
	@test isa(s, StatsBase.SummaryStats)
	@test s.min    ≈ 0.0
	@test s.max    ≈ 1.0
	@test s.mean   ≈ 0.9286052603644784
	@test s.median ≈ 0.9963369963369964
	@test s.q25    ≈ 0.9725274725274725
	@test s.q75    ≈ 0.9990842490842491

finally
	seqClose(f)
end




## Test: the calculation of covariance matrix for PCA

f = seqOpen(seqExample(:kg))
println("Calculation of covariance matrix for PCA")

try
	# initialize the covariance matrix
	cov1 = Any[0.0]
	# apply a user-defined function
	seqApply(f, "#dosage", cov1) do geno::Matrix{UInt8}, cov1::Vector{Any}
		# calculate allele frequencies
		N = size(geno, 2); af = Vector{Float64}(N)
		for i in 1:N
			g = geno[:,i]; af[i] = mean(g[g .!= 0xFF])
		end
		af *= 0.5
		# normalized by allele frequencies
		g = Matrix{Float64}(size(geno))
		for i in 1:N
			g[:,i] = (geno[:,i] - 2*af[i]) / sqrt(af[i] * (1 - af[i]))
		end
		# correct missing genotypes
		g[isnan(g)] = 0.0; g[geno .== 0xFF] = 0.0
		# update the cov matrix
		cov1[1] += g * g'
	end
	# scale the matrix
	cov1 = cov1[1] * size(cov1[1], 1) / sum(diag(cov1[1]))

	s = summarystats([ cov1... ])
	@test isa(s, StatsBase.SummaryStats)
	@test s.min    ≈ -0.0667281080577368
	@test s.max    ≈ 3.3600075014863386
	@test_approx_eq_eps s.mean 0.0 1e-12
	@test s.median ≈ -0.006934683064608682
	@test s.q25    ≈ -0.030999663540459556
	@test s.q75    ≈ 0.016731842043406

finally
	seqClose(f)
end


