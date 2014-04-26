module Dets

export det

using ..Common.Generator


code_det_kernel(A::Matrix) = code_det_kernel(A, "subdet")
function code_det_kernel(A::Matrix, sdname::String)
    n = size(A,1)
    @assert n == size(A,2)

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
        if j < 2; continue; end

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

    code, subdets[end]
end

# ---- det(rows...) ----
for NC = 1:8
    generators = [symbol("v$k") for k in 1:NC]
#    args = [:($generator::Generator{$NC}) for generator in generators]
    args = [:($generator) for generator in generators]

    A  = [ :($generator[$i]) for i=1:NC, generator in generators ]
    code, d = code_det_kernel(A, "d")    

    @eval function det($(args...))
        $(code...)
        $d        
    end
end

end # module
