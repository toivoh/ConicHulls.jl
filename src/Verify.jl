module Verify

export create_nb_dict, get_nb

using ..Common


immutable Link{G}
    generators::Vector{G}
    even::Bool

    function Link(generators::Vector, even::Bool, k::Int)
        new(except_index(generators, k), even $ isodd(length(generators)-k))
    end
end
Base.isequal(l1::Link, l2::Link) = (l1.generators == l2.generators && l1.even == l2.even)
Base.hash(link::Link) = hash(link.generators) $ link.even

function add_links!{F<:Facet,L<:Link}(links::Dict{L,F}, facets)
    for facet in facets
        generators, even = get_canonical_winding(facet)
        for k in 1:nfacet(facet)
            link = L(generators, even, k)
            if haskey(links, link); error("Link already exists!"); end
            links[link] = facet
        end
    end
end

function create_nb_dict{F<:Facet}(::Type{F}, facets)
    links = Dict{Link{gtype(F)}, F}()
    add_links!(links, facets)
    links
end

function get_nb{F<:Facet,L<:Link}(links::Dict{L,F}, facet::F, g)
    generators, even = get_canonical_winding(facet)
    k = indexof(generators, g)
    link = L(generators, !even, k)
    links[link]
end


end # module
