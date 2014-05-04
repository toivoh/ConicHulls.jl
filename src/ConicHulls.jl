module ConicHulls

export Facet, nconic, nface, ftype, gtype
export ConicHull, create_hull, add!, verify

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

Hulls.create_hull(NC::Int, generators) = create_hull(hulltype(NC), generators)


end # module
