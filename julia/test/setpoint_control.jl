function c2d(Ac::AbstractArray, Bc::AbstractArray, h::Float64)
    A = exp(Ac*h)
    B = Ac\(A-I)*Bc
    return A,B
end      

# linear harmonic oscillator 
k = 1.0     # spring constant 
m = 0.05     # mass
b = 0.00     # damping coefficient 
Ac = [0 1; -k/m -b/m]
Bc = [0; 1/m]   
 
h = 0.01    # discretization 
A,B = c2d(Ac,Bc,h)      

function linear_harmonic_oscillator(x,u)
    return A*x+B*u
end  

function rk4(x, u, h, dynamics)
    k1 = dynamics(x, u)
    k2 = dynamics(x + 0.5h*k1, u)
    k3 = dynamics(x + 0.5h*k2, u)
    k4 = dynamics(x + h*k3, u)
    return x + h/6*(k1 + 2k2 + 2k3 + k4)
end  

Nsim = 2000
thist=0:h:(Nsim-1)*h
xhist = zeros(2,Nsim)
xhist[:,1] = [2π;0]
for i in 1:Nsim-1
    xhist[:,i+1] = linear_harmonic_oscillator(xhist[:,i],0.0)
end 
       
fig = Figure()
ax = Axis(fig[1,1])
lines!(ax,thist,xhist[1,:])
fig    
       
lho = LTImodel(A,B) 

begin 
    Q = [1;0] 
    R = [0.5] 
    xf = [0;0] 
    hx = 80

    Δumax= [Inf] 
    Δumin= [Inf] 
    umax = [10] 
    umin = [-10] 
    xmax = [1.0;20] 
    xmin = [-100;-Inf] 

    c = Constraints(lho,Δumax,Δumin,umax,umin,xmax,xmin,hx);
    
    ctrl = MPCincremental(lho,Q,Q,R,xf,hx,c);   

    Nsim = 500
    thist=0:h:(Nsim-1)*h
    xhist = zeros(2,Nsim)
    xhist[:,1] = [-2π;0]  
    uhist = zeros(1,Nsim-1)
    for i in 1:Nsim-1
        u = get_control(xhist[:,i], ctrl)
        uhist[:,i]=u
        xhist[:,i+1] = linear_harmonic_oscillator(xhist[:,i],u[1])
    end    
    
    begin 
        fig = Figure()
        ax1 = Axis(fig[1,1])
        lines!(ax1,thist[1:end-1],uhist[1,:]) 
        ax2 = Axis(fig[2,1])
        lines!(ax2,thist[1:end],xhist[1,:]) 
        ax3 = Axis(fig[3,1])
        lines!(ax3,thist[1:end],xhist[2,:]) 
        fig       
    end  
end     

      
 