module Common

export Generator, Facet, nconic, nfacet
export except_index, indexof


abstract Generator{NC}
abstract Facet{NF}

nconic{NC}(::Type{Generator{NC}}) = NC
nconic{NC}(::Generator{NC})       = NC
nconic{NF}(::Type{Facet{NF}})     = NF + 1
nconic{NF}(::Facet{NF})           = NF + 1

nfacet(x) = nconic(x) - 1

except_index{T}(x::Vector{T}, k::Int) = [x[1:(k-1)]... x[(k+1):end]...]

function indexof(v::Vector, x)
    for (k, y) in enumerate(v); if y == x;  return k;  end; end
    error("indexof: x = $x not found in v=$v")
end

end # module
