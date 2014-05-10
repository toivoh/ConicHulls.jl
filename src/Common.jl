module Common

export Generator, Facet, nconic, nfacet, gtype
export get_canonical_winding
export except_index, indexof, isevenperm, isoddperm


abstract Generator{NC}
abstract Facet{NF}

nconic{NC}(::Type{Generator{NC}}) = NC
nconic{NC}(::Generator{NC})       = NC
nconic{NF}(::Facet{NF})           = NF + 1

nfacet(x) = nconic(x) - 1

gtype(::Facet)                 = error("Unimplemented")
get_canonical_winding(::Facet) = error("Unimplemented")

except_index{T}(x::Vector{T}, k::Int) = [x[1:(k-1)]... x[(k+1):end]...]

function indexof(v::Vector, x)
    for (k, y) in enumerate(v); if y == x;  return k;  end; end
    error("indexof: x = $x not found in v=$v")
end

isoddperm(p) = !isevenperm(p)
function isevenperm{T}(p::Vector{T})
    @assert isperm(p)
    n = length(p)
    used = falses(n)
    even = true
    for k = 1:n
        if used[k]; continue; end
        # Each even cycle flips even (an odd number of times)
        used[k] = true
        j = p[k]
        while !used[j]
            used[j] = true
            j = p[j]
            even = !even
        end
    end
    even
end

end # module
