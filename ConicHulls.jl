module ConicHulls

include("Common.jl")
include("Dets.jl")
include("Primitives.jl")
include("Hulls.jl")

using .Primitives
using .Hulls


hulltype(NC) = Hulls.hulltype(AGen{NC,Int})

function create_simplex_hull(NC::Int)
    Hulls.create_simplex_hull(hulltype(NC))
end


end # module
