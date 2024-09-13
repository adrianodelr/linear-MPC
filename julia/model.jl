mutable struct LTImodel
    A::AbstractArray
    B::AbstractArray
    n::Int
    m::Int 
    x0::AbstractArray
    u0::AbstractArray
    function LTImodel(A::AbstractArray,B::AbstractArray)
        n = size(A,1)
        m = size(B,2)
        return new(A,B,n,m,zeros(n),zeros(m))
    end 
end 


