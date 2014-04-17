module ConicHulls

export hulltype, create_simplex_hull

include("Common.jl")
include("Dets.jl")
include("Primitives.jl")
include("Hulls.jl")

using .Primitives
using .Hulls


Hulls.hulltype(NC) = hulltype(AGen{NC,Int})

function Hulls.create_simplex_hull(NC::Int)
    create_simplex_hull(hulltype(NC))
end


end # module
