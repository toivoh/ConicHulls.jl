module ConicHulls

export Facet, nconic, nface, ftype, gtype
export hulltype, create_simplex_hull, add!, verify

include("Common.jl")
include("Dets.jl")
include("Primitives.jl")
include("Verify.jl")
include("Hulls.jl")
include("RefHull.jl")

using .Common
using .Primitives
using .Hulls


Hulls.hulltype(NC) = hulltype(AGen{NC,Int})

function Hulls.create_simplex_hull(NC::Int)
    create_simplex_hull(hulltype(NC))
end


end # module