function [K,kprime,marginal_k, StDist, Gamma] = K_Agg(R,w,par,mpar,P,meshes,gri)
    % This function calculates the aggregate supply of funds for a given
    % interest rate R and wage rate w. 
    if par.gamma ==1
        mutil    = @(c) 1./c;  % Marginal utility
        invmutil = @(mu) 1./mu;% inverse marginal utility
    else
        mutil    = @(c) 1./(c.^par.gamma); % Marginal utility
        invmutil = @(mu) 1./(mu.^(1./par.gamma)); % inverse marginal utility
    end
    par.r    = R;
    meshes.z = meshes.z*w; % take into account wages
    Cold     = (meshes.z  + par.r*meshes.k); %Initial guess for consumption policy: roll over assets
    Cold(Cold<0.01) = 0.01; % to avoid negative consumption
    
    distEG   = 1; % Initialize Distance
    while distEG>mpar.crit
        C      = EGM(Cold,mutil,invmutil,par,mpar,P,meshes,gri); % Update consumption policy by EGM
        distEG = max(abs(C(:)-Cold(:))); % Calculate Distance
        Cold   = C; % Replace old policy
        Cold(Cold<0.01) = 0.01; % to avoid negative consumption
    end
    [~,kprime] = EGM(C, mutil, invmutil,par,mpar,P,meshes,gri);
    
    [Gamma, StDist] = Young(kprime,gri,mpar,P);
    marginal_k      = sum(reshape(StDist,[mpar.nk, mpar.nz]),2);
    K               = dot(marginal_k,gri.k);
    end