## Solve the Consumption Savings Model using Value Function Iteration
using Plots
include("Functions/VFI_update_spline.jl")

## 1. Define parameters

# Numerical parameters
# `nk`: Number of points on the asset grid
# `nz`: Number of points on the log-productivity grid
# `crit`: Numerical precision
# `maxk`: Maximum assets
# `mink`: Minimum assets (equal to borrowing limit)
mpar = (nk = 100, nz = 2, crit = 1.0e-5, maxk = 6, mink = -9/4)
println("Numerical parameters")
println(mpar) # Display numerical parameters

# Economic Parameters
# `r`: Real Rate
# `γ`: Coefficient of relative risk aversion
# `β`: Discount factor
# `b`: borrowing limit
par = (r = 4/90, γ = 1.0, β = 0.95, b = mpar.mink)
println("Economic parameters")
println(par) # Display economic parameters

## 2. Generate grids, Meshes and Income
#Define asset grid on log-linearspaced
gri   = (
            k = xxx,
            z = [1/9; 10/9]
        )
Π     = [3/5 2/5; 4/90 86/90] # Transition matrix for income


## 3. Define utility functions

if par.γ == 1.0
    util(c)     = xxx # Utility
    mutil(c)    = xxx  # Marginal utility
    invmutil(mu) = xxx # inverse marginal utility
else
    util(c)      = xxx # Utility
    mutil(c)     = xxx # Marginal utility
    invmutil(mu) = xxx # inverse marginal utility
end


## 4b. Value Function Iteration (on-grid)
# Meshes of capital and productivity
meshes = (
            kprime = [kP for kP in gri.k, k in gri.k, z in gri.z],
            k = [k for kP in gri.k, k in gri.k, z in gri.z],
            z = [z for kP in gri.k, k in gri.k, z in gri.z]
        )
Y = meshes.z + meshes.k .*(1+par.r) # Cash at hand (Labor income plus assets cum dividend)
timer = time() # Start timer
V               = zeros(xxx,mpar.nz) # Initialize Value Function
distVF_on       = ones(1) # Initialize Distance
iterVF          = 1 # Initialize Iteration count
C_hat           = xxx
U               = fill(-1.0e10,size(C_hat))
U[C_hat .> 0]   = xxx(C_hat[C_hat .> 0])
while distVF_on[iterVF] > mpar.crit # Value Function iteration loop: until distance is smaller than crit.
    # Update Value Function using on-grid search
    local EV              =  xxx         
    local V_hat           = U + repeat(reshape(par.β*EV,mpar.nk,1,mpar.nz), outer=(1,mpar.nk,1))
    local pol_ind         = argmax(xxx,dims=1) # Optimize given cont' value
    local Vaux            = xxx
    local Vnew            = dropdims(Vaux,dims=1)
    global kprime_on       = [gri.k[cix[1]] for cix in dropdims(pol_ind,dims=1)]
    local dd              = maximum(abs.(Vnew-V)) # Calculate distance between old guess and update
    global V               = Vnew # Update Value Function
    global iterVF          += 1 # Count iterations
    append!(distVF_on,dd)   # Save distance
end
V_on       = reshape(V,mpar.nk,mpar.nz)
time1 = time() - timer # Save Time used for VFI


## 4b. Value Function Iteration (off-grid)
# Meshes of capital and productivity
meshes = (
            k = [k for k in gri.k, z in gri.z],
            z = [z for k in gri.k, z in gri.z]
         )
Y = meshes.z + meshes.k*(1+par.r) # Cash at hand (Labor income plus assets cum dividend)
timer = time() # Start timer
V    = zeros(mpar.nk,mpar.nz) # Initialize Value Function
distVF = ones(1) # Initialize Distance
iterVF = 1 # Initialize Iteration count
while distVF[iterVF]>mpar.crit # Value Function iteration loop: until distance is smaller than crit.
    # Update Value Function using off-grid search
    global Vnew,kprime = VFI_update_spline(xxx,xxx,xx,par,mpar,gri,Π) # Optimize given cont' value
    local dd          = maximum(abs.(Vnew-V)) # Calculate distance between old guess and update

    global V          = Vnew # Update Value Function
    global iterVF     += 1 # Count iterations
    append!(distVF,dd)   # Save distance
end
V      = reshape(V,mpar.nk,mpar.nz)
time2 = time() - timer # Save Time used for VFI

## 5. Plot Value and Policy Functions from VFI

figure1 = plot(gri.k,V_on,linestyle = :dash,labels=["low productivity (on)" "high productivity (on)"],legend = :bottomright) # Plot value functions
plot!(gri.k,V,labels=["low productivity (off)" "high productivity (off)"]) # Plot value functions
title!("Value Function from VFI") # Title and legend of the graph
xlabel!("assets")

figure2 = plot(gri.k,kprime_on,linestyle = :dash,labels=["low productivity (on)" "high productivity (on)"],legend = :topleft) # Plot policy functions
plot!(gri.k,kprime,labels=["low productivity (off)" "high productivity (off)"]) # Plot policy functions
plot!(gri.k,gri.k,linestyle = :dot,linewidth = 2,linecolor= :black,labels = nothing) # Add 45° line
title!("Policy Function from VFI") # Title of the graph
xlabel!("assets")
ylabel!("saving")

## 6. Display time and convergence stats
println("Time for solution")
println("on grid:",time1)
println("off grid: ",time2)

figure3 = plot(1:length(distVF_on)-1,distVF_on[2:end],yaxis = :log,label="on grid")
plot!(1:length(distVF)-1,distVF[2:end],label="off grid")
title!("distance after iterations")