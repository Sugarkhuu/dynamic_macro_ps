using CubicSplines, Optim

"""
    VFI_update_spline(V,Y,util,par,mpar,gri,Π)

Update the value function (one VFI iteration) for the consumption-savings problem.

    * `V` (dimensions: k x z) is the old value function guess.
    * `Y` (dimensions: k x z) is a matrix of cash at hand
    * `util` is the felicity function.
    * `par` and `mpar` are structures containing economic and numerical  parameters.
    * `gri` is the asset and productivity grid.
    * `Π` (dimensions: z x z') is the transition probability matrix.

Returns

    * `Vnew`: (k x z) optimal value function
    * `kprime`: (k x z) asset policy
"""
function VFI_update_spline(V,Y,util,par,mpar,gri,Π)
    V      = reshape(V,mpar.nk,mpar.nz) # make sure that V has the right format dim1: k, dim2:z
    kprime = zeros(size(V)) # allocate policy matrix
    Vnew   = zeros(size(V)) # allocate new value matrix
    EV     =  xxx   # Calculate expected continuation value

    for zz = 1:mpar.nz # loop over Incomes
        ev_int = CubicSpline(xxx,xxx) # Prepare interpolant
        for kk=1:mpar.nk # loop of Assets
            f(kpri)       = -util(xxx)-par.β*ev_int[xxx] # Define function to be minimized
            optimres        = optimize(f,par.b,min(xxx,xxx)) # Find minimum of f for savings between par.b and min(Y[kk,zz],maxk)
            Vnew[kk,zz]   = - Optim.minimum(optimres) # Save value
            kprime[kk,zz] = Optim.minimizer(optimres) # Save policy
        end
    end
    return Vnew, kprime
end