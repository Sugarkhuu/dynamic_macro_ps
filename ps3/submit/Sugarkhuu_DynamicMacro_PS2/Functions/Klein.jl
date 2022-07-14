using LinearAlgebra

"""
    Klein(A,B,mpar)

Compute state-to-control mapping and state transition
from first-order system dynamics

    * `A`: Jacobian of system wrt. future variables
    * `B`: Jacobian of system wrt. present variables
    * `mpar`: parameter structure

Returns

    * `S2C`: state-to-control mapping
    * `S2Sprime`: state transition
"""
function Klein(A,B,mpar)
    Haha = schur(complex(A),-complex(B))
    
    relev = abs.(diag(Haha.S))./abs.(diag(Haha.T))
    ll    = sort(relev)
    slt   = relev .>= 1.0
    nk    = sum(slt)
    
    if nk>mpar.nstates
        @warn "The Equilibrium is Locally Indeterminate!"
    elseif nk<mpar.nstates
        @warn "No Local Equilibrium Exists!"
    end
    ordschur!(Haha,slt)
    
    z21=Haha.Z[nk+1:end,1:nk]
    z11=Haha.Z[1:nk,1:nk]
    s11=Haha.S[1:nk,1:nk]
    t11=Haha.T[1:nk,1:nk]
    
    # Checks
    if rank(z11)<nk
        @warn "invertibility condition violated"
    end
    z11i=z11\I(nk)
    S2C=real(z21*z11i) #States2Controls
    S2Sprime=real(z11*(s11\t11)*z11i) #LOM states
    return S2C,S2Sprime
end