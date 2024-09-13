struct QPweights
    Q::AbstractArray
    R::AbstractArray
    function QPweights(Q::Vector, Qf::Vector, R::Vector, hx::Int)
        Qf = diagm(Qf)
        n = length(Q)
        Q = diagm(repeat(Q, inner=[1,1], outer=[1,hx]) |> vec)
        R = diagm(repeat(R, inner=[1,1], outer=[1,hx]) |> vec)
        Q[end-n+1:end,end-n+1:end] .= Qf    
        return new(Q,R)
    end 
end 

struct PredictionMatrices
    Ā::AbstractArray
    B̄::AbstractArray
    function PredictionMatrices(Aa::AbstractArray, Ba::AbstractArray, hx::Int)
        n = size(Aa,1)
        m = size(Ba,2)        
        Ā=zeros(hx*n,n);
        B̄=zeros(hx*n,hx*m);
        for i=1:hx
            Ā[(1+(i-1)*n):i*n,1:n]=Aa^i;
            for j=1:hx            
                if(i>=j)                
                    B̄[(1+(i-1)*n):i*n,(1+(j-1)*m):j*m] = Aa^(i-j)*Ba;                
                end
            end
        end
        return new(Ā,B̄)
    end 
end 

struct Constraints 
    Δumax::AbstractArray
    Δumin::AbstractArray
    umax::AbstractArray
    umin::AbstractArray
    xmax::AbstractArray
    xmin::AbstractArray
    LΔ::AbstractArray
    pMat::PredictionMatrices
    function Constraints(model::LTImodel,Δumax::AbstractArray,Δumin::AbstractArray,umax::AbstractArray,umin::AbstractArray,xmax::AbstractArray,xmin::AbstractArray,hx::Int)
        m = length(umin)
        Δumax = repeat(Δumax,hx)
        Δumin = repeat(Δumin,hx)
        umax  = repeat(umax,hx)
        umin  = repeat(umin,hx)
        xmax  = repeat(xmax,hx)
        xmin  = repeat(xmin,hx)
        LΔ    = kron(ones(hx,hx) |> tril ,I(m))
        pMat  = PredictionMatrices(model.A,model.B,hx)
        return new(Δumax,Δumin,umax,umin,xmax,xmin,LΔ,pMat)
    end 
end

function build_QP(x::Vector, xref::Vector, w::QPweights, m::PredictionMatrices)
    P = m.B̄'*w.Q*m.B̄+w.R
    q = m.B̄'*w.Q * (m.Ā*x - xref)
    return P,q
end 

mutable struct MPCincremental
    model::LTImodel
    predMat::PredictionMatrices
    hx::Int
    na::Int 
    QPw::QPweights
    xref::Vector
    c::Constraints 
    solver 
    function MPCincremental(model::LTImodel,Q::Vector,Qf::Vector,R::Vector,xgoal::Vector,hx::Int, c::Constraints)
        n,m = model.n, model.m
        Aa = [model.A model.B; zeros(m,n) I(m)]
        Ba = [model.B; I(m)]
        na = n+m
        QPw = QPweights([Q;zeros(m)], [Qf;zeros(m)], R, hx)
        pMat = PredictionMatrices(Aa, Ba, hx)
        solver = OSQP.Model()
        xref = repeat([xgoal;zeros(m)],hx)
        lb,ub,A = build_constraints(zeros(n),c, model, hx)
        P,q = build_QP(zeros(na), xref, QPw, pMat)
        OSQP.setup!(solver, P=sparse(P), q=vec(q), A=SparseMatrixCSC(A),l=lb, u=ub, verbose=false,warm_start=true)  
        return new(model, pMat, hx, na, QPw, xref, c, solver)
    end 
end 

function get_control(x::Vector, ctrl::MPCincremental)
    P,q = build_QP([x;ctrl.model.u0], ctrl.xref, ctrl.QPw, ctrl.predMat)
    lb,ub,A = build_constraints(x,ctrl.c, ctrl.model, ctrl.hx)
    OSQP.update!(ctrl.solver, Px=triu(sparse(P)).nzval, q=vec(q), l=lb, u=ub, Ax=sparse(A).nzval)
    Δu = OSQP.solve!(ctrl.solver).x
    u = ctrl.model.u0 + Δu[1:ctrl.model.m]
    ctrl.model.u0 = u 
    ctrl.model.x0 = x
    return u 
end 

function build_constraints(x::Vector, c::Constraints, m::LTImodel, hx::Int)
    Lu = repeat(m.u0,hx)
    lb = [c.umin-Lu;c.xmin-c.pMat.Ā*x-c.pMat.B̄*Lu]
    ub = [c.umax-Lu;c.xmax-c.pMat.Ā*x-c.pMat.B̄*Lu]
    A  = [c.LΔ;c.pMat.B̄*c.LΔ]
    return lb,ub,A
end 
