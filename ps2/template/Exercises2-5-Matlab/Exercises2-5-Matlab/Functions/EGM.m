%% Policy update by EGM
function [C,Kprime] = EGM(C,mutil,invmutil,par,mpar,P,meshes,gri)
    %% This function iterates forward the consumption policies for the consumption Savings
    % model using the EGM method. C (k x z) is the consumption policy guess. MUTIL and INVMUTIL are
    % the marginal utility functions and its inverse. PAR and MPAR are parameter structures.
    % P is the transition probability matrix. MESHES and GRI are meshes and grids for income
    % (z) and assets (k).

    mu     = xxx; % Calculate marginal utility from c'
    emu    = xxx;     % Calculate expected marginal utility
    Cstar  = xxx;     % Calculate cstar(m',z)
    Kstar  = xxx; % Calculate mstar(m',z)
    Kprime = xxx; % initialze Capital Policy

    for z=1:mpar.nz % For massive problems, this can be done in parallel
        % generate savings function k(z,kstar(k',z))=k'
        Savings     = griddedInterpolant(xxx,xxx,'spline');
        Kprime(:,z) = Savings(xxx);     % Obtain k'(z,k) by interpolation
        BC          = xxx < Kstar(1,z); % Check Borrowing Constraint
        % Replace Savings for HH saving at BC
        Kprime(BC,z)= par.b; % Households with the BC flag choose borrowing contraint
    end
    % generate consumption function c(z,k^*(z,k'))
    C          = xxx; %Consumption update
end
