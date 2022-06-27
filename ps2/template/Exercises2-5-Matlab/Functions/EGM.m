%% Policy update by EGM
function [C,Kprime] = EGM(C,mutil,invmutil,par,mpar,P,meshes,gri)
    %% This function iterates forward the consumption policies for the consumption Savings
    % model using the EGM method. C (k x z) is the consumption policy guess. MUTIL and INVMUTIL are
    % the marginal utility functions and its inverse. PAR and MPAR are parameter structures.
    % P is the transition probability matrix. MESHES and GRI are meshes and grids for income
    % (z) and assets (k).

    mu     = mutil(C); % Calculate marginal utility from c'
    emu    = mu*P';     % Calculate expected marginal utility
    Cstar  = invmutil((1+par.r)*par.beta*emu);     % Calculate cstar(m',z)
    Kstar  = (Cstar + meshes.k - meshes.z)/(1+par.r); % Calculate mstar(m',z)
    Kprime = meshes.k; % initialze Capital Policy
keyboard
    for z=1:mpar.nz % For massive problems, this can be done in parallel
        % generate savings function k(z,kstar(k',z))=k'
        Savings     = griddedInterpolant(Kstar(:,z),gri.k,'spline');
        Kprime(:,z) = Savings(gri.k);     % Obtain k'(z,k) by interpolation
        BC          = gri.k < Kstar(1,z); % Check Borrowing Constraint
        % Replace Savings for HH saving at BC
        Kprime(BC,z)= par.b; % Households with the BC flag choose borrowing contraint
    end
    disp("he")
    % generate consumption function c(z,k^*(z,k'))
    C          = (1+par.r)*meshes.k+meshes.z-Kprime; %Consumption update
end
