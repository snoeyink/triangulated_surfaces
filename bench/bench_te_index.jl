using BenchmarkTools

#include("src/TriangulatedSurfaces.jl")
#using .TriangulatedSurfaces

function test_int(inttype, n::Int)
	# counter: length-n vector of the requested integer type, zero-initialized
	counter = zeros(inttype, n)
	# digits: small fixed buffer of 20 entries of the same integer type
	digits = zeros(inttype, 20)

    i = n
    while i > 0
        # increment the counter vector as if it were a number in base 20
        @inbounds while i > 0
            if counter[i] > 0
                digits[counter[i]] -= 1
            end
            if counter[i] >= 20
                counter[i] = 1 
                digits[1] += 1
                i -= 1
            else
                counter[i] += 1
                digits[counter[i]] += 1
                i = n
                break
            end
        end

    end
	return counter, digits
end

function benchmark_te()
    for t in (Int64, Int32, UInt16, UInt8, UInt64, UInt32, UInt16, UInt8)
        println()
        println(t)
        @btime begin c,d = test_int($t, 5); sum(d) end
    end
    return nothing
end


