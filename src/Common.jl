module Common

export Generator, Facet, nconic, nfacet, gtype
export get_canonical_winding
export except_index, indexof, isevenperm, isoddperm

export dominates, antidominates
export Face, facesof, nbof, opposite, set_opposite!, replace_link!, set_nb!

abstract Generator
abstract Facet

nconic(::None) = error("Unimplemented")

nfacet(x) = nconic(x) - 1

dot(f::None, g::None) = error("Unimplemented")
dominates(    g::Generator, f::Facet) = dot(f, g) > 0
antidominates(g::Generator, f::Facet) = dot(f, g) < 0


gtype(::None)                 = error("Unimplemented")
get_canonical_winding(::None) = error("Unimplemented")

except_index{T}(x::Vector{T}, k::Int) = [x[1:(k-1)]..., x[(k+1):end]...]

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


set_opposite!(f::None, new_link::None, g::None)    = error("Unimplemented")
replace_link!(f::None, new_link::None, link::None) = error("Unimplemented")

# ---------------------------------- Face ------------------------------------

immutable Face{F<:Facet}
    parent::F
    k::Int
end
Face{F<:Facet}(parent::F, nb::F) = Face{F}(parent, indexof(parent, nb))
Face{F<:Facet}(parent::F, g::Generator) = Face{F}(parent, indexof(parent, g))

parentof(f::Face) = f.parent
nbof(f::Face)     = f.parent.links[f.k]
opposite(f::Face) = f.parent.generators[f.k]
flip(f::Face)     = (nb = nbof(f); Face(nb, indexof(nb, f)))

set_nb!{F}(f::Face{F}, newnb::F) = (f.parent.links[f.k] = newnb)


# ---------------------------------- Faces -----------------------------------

immutable Faces{F<:Facet}; parent::F; end

Base.length(faces::Faces) = nfacet(faces.parent)
Base.start(faces::Faces) = 1
Base.next{F}(faces::Faces{F}, k) = (Face{F}(faces.parent, k), k+1)
Base.done(faces::Faces, k) = k == (nfacet(faces.parent)+1)
Base.eltype{F}(faces::Faces{F}) = Face{F}

facesof(f::Facet) = Faces(f)


end # module
