module RefHull
export verify_hull

using ConicHulls.Dets.incomplete_det

function eliminate_duplicates{T<:Integer}(generators::Vector{Vector{T}})
    found = Set{NTuple{length(generators[1]), T}}()
    gs = Vector{T}[]
    for g in generators
        gd = div(g, gcd(g...))
        key = tuple(gd...)
        if !(key in found)
            push!(gs, g)
            push!(found, key)
        end
    end
    gs
end

type SortedFacet{T}
    generators::Vector{Vector{T}}
    positive::Bool
    plane::Vector{T}
end

function enumerate_planes{T}(NC::Int, generators::Vector{Vector{T}})
    NF = NC - 1
    extreme = ObjectIdDict()
    facets = SortedFacet{T}[]
    for ks in combinations(1:length(generators), NF)
        gs = generators[ks]
        plane = incomplete_det(gs...)

        if all(plane .== 0); continue; end
    
        positive = negative = npos = nneg = 0
        for (k,g) in enumerate(generators)
            if k in ks; continue; end
            
            d = dot(plane, g)
            if d > 0;     positive = k; npos += 1;                 
            elseif d < 0; negative = k; nneg += 1;
            end
            if (npos > 1) && (nneg > 1); break; end
        end
        if npos == 1; extreme[generators[positive]] = true; end
        if nneg == 1; extreme[generators[negative]] = true; end
        if npos == nneg == 0
            error("Generators do not span the the space")
        end
        if npos == 0 || nneg == 0
            positive = (npos == 0)
            push!(facets, SortedFacet(gs, positive, positive ? plane : -plane))
        end
    end
    # Filter out non-extreme generators
    gs = Vector{T}[]
    for g in generators
        if haskey(extreme, g); push!(gs, g); end
    end
    # Filter out facets with non-extreme generators
    fs = SortedFacet{T}[]
    for facet in facets
        if all([haskey(extreme, g) for g in facet.generators])
            push!(fs, facet)
        end
    end
    gs, fs
end

function filter_coplanar_facets{T}(generators::Vector{Vector{T}}, facets::Vector{SortedFacet{T}})
    NC = length(generators[1])
    order = ObjectIdDict()
    for (k,g) in enumerate(generators); order[g] = k; end

    fs = SortedFacet{T}[]
    for facet in facets
        j = first(find(facet.plane .!= 0))
        dominated = false
        for g in generators
            if g in facet.generators; continue; end
            
            d = dot(facet.plane,g)
            @assert d <= 0
            if d < 0; continue; end

            # Coplanar
            gs = Vector{T}[facet.generators..., g]

            ks = Int[order[g] for g in facet.generators]
            push!(ks, order[g])
            
            for i in sortperm(ks, rev=true)
                p = incomplete_det(gs[1:(i-1)]...,gs[(i+1):end]...)
                if all(p .== 0); continue; end

                sameway = (p[j] > 0) == (facet.plane[j] > 0)
                if (sameway == iseven(NC-i)) == facet.positive
                    dominated = true; 
                else
                    @assert i < NC # i == NC ==> should always be dominated
                end
                break
            end
            if dominated; break; end
        end
        if !dominated; push!(fs, facet); end
    end
    fs
end

function conichull{T<:Integer}(NC::Int, generators::Vector{Vector{T}})
    # Verify input data
    for g in generators; @assert length(g) == NC; end

    generators         = eliminate_duplicates(generators)
    generators, facets = enumerate_planes(NC, generators)
    facets             = filter_coplanar_facets(generators, facets)

    generators, facets
end

using ConicHulls: Common, Hulls

function verify_hull(hull::ConicHull, generators::Vector) 
    orig_generators = generators
    validate(hull)

    gmap = ObjectIdDict()
    for g in generators; gmap[g.x] = g; end
    generators = [g.x for g in generators]
    generators, facets = RefHull.conichull(nconic(hull), generators)

    hfacets = {}
    for facet in hull.facets
        gs, even = get_canonical_winding(facet)
        push!(hfacets, (gs, even))
    end
    sort!(hfacets, 
          lt=(w1,w2)->(
              for (g1,g2) in zip(w1[1], w2[1]);
                  if g1.id < g2.id; return true; end;
                  if g2.id > g2.id; return false; end;
              end;
              return false;))

    try
        @assert length(hull.facets) == length(facets)
        facetset = Set()
        for facet in facets
            push!(facetset, tuple([gmap[g] for g in facet.generators]..., facet.positive))
        end
        for facet in hull.facets
            gs, even = get_canonical_winding(facet)
            key = tuple(gs..., even)
            @assert key in facetset
        end
    catch e
        @show length(hull.facets), length(facets)

        println("hull.facets:")
#        for facet in hull.facets; println("  ", facet.generators); end
        for (gens, pos) in hfacets; println("  ", gens, ": ", pos); end
        println("facets:")
        for facet in facets; println("  ", [gmap[g] for g in facet.generators], ": ", facet.positive); end
        println("generators:")
        for g in orig_generators; println("  ", g, ": ", g.x); end
        println("filtered generators:")
        for g in generators; println("  ", gmap[g]); end

        rethrow(e)
    end
end

end # module
