module Dets

export det

using ..Common.Generator


code_det_kernel(A::Matrix) = code_det_kernel(A, "subdet")
function code_det_kernel(A::Matrix, sdname::String)
    n, m = size(A)
    @assert m <= n
    @assert m >= 1

    code = {}
    N = 1<<n
#    subdets = {gensym("subdet$k") for k=1:(N-1)}
#    subdets = {symbol("$sdname$(bin(k))") for k=1:(N-1)}
    subdets = {symbol(string(sdname, 
                             [(k&(1<<(j-1))!=0 ? j : "") for j=1:n ]...))
               for k=1:(N-1)}

    for i=1:n;   subdets[1<<(i-1)] = A[i,1];  end    
    for mask = 1:(N-1)
        j = count_ones(mask)
        if j < 2 || j > m; continue; end

        terms = {}
        negate = false
        for i in n:-1:1
            bit = 1<<(i-1)
            if (mask & bit) == 0; continue; end
            
            term = A[i,j]
            if j > 1
                subdet = subdets[mask & ~bit]
                term = :($subdet*$term)
            end
            if negate;  term = :(-$term);  end
            push!(terms, term)
            negate = !negate
        end
        push!(code, :($(subdets[mask]) = +($(terms...))))
    end    

    code, subdets
end

function incomplete_det_components(n::Int, subdets::Vector)
    ds = {}
    for k=1:n
        d = subdets[(1<<n)-1-(1<<(k-1))]
        if isodd(n-k);  d = :(-$d); end
        push!(ds, d)
    end
    ds
end


incomplete_det() = [1]

# ---- det(rows...), incomplete_det(rows...) ----
for NC = 1:8
    generators = [symbol("v$k") for k in 1:NC]
#    args = [:($generator::Generator{$NC}) for generator in generators]
    args = [:($generator) for generator in generators]

    A  = [ :($generator[$i]) for i=1:NC, generator in generators ]
    code, ds = code_det_kernel(A, "d")
    d = ds[end]

    @eval function det($(args...))
        $(code...)
        $d        
    end

    if NC > 1
        code, ds = code_det_kernel(A[:,1:(end-1)], "d")
        ds = incomplete_det_components(NC, ds)
        iargs = args[1:(end-1)]
        @eval function incomplete_det($(iargs...))
            $(code...)
            [$(ds...)]
        end
    end
end

end # module
