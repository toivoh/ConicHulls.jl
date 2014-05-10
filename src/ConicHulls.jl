module ConicHulls

export Facet, nconic, nface, ftype, gtype
export ConicHull, create_hull, init_hull!, add!, verify

include("Common.jl")
include("Dets.jl")
include("Primitives.jl")
include("Verify.jl")
include("Hulls.jl")
include("RefHull.jl")

using .Common
using .Primitives
using .Hulls


Hulls.hulltype(NC) = Hulls.hulltype(AGen{Int})

Hulls.create_hull(NC::Int) = create_hull(NC, Hulls.hulltype(NC))
Hulls.create_hull(NC::Int, generators) = create_hull(NC, Hulls.hulltype(NC), generators)


end # module
