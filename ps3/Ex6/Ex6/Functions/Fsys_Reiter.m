function F = Fsys_Reiter(Xhat, XhatPrime, XSS, mpar, par, id,gri,meshes,P)
% expand X around steady state

X      = exp(Xhat + XSS);
XPrime = exp(XhatPrime+ XSS);

Dist            = zeros(mpar.nk, mpar.nz);
Dist(1:end-1)   = Xhat(id.Dev.Dist) + XSS(id.SS.Dist);
Dist(end)       = 1 - sum(Dist(:));
Dist            = Dist(:);

DistPrime            = zeros(mpar.nk, mpar.nz);
DistPrime(1:end-1)   = XhatPrime(id.Dev.Dist) + XSS(id.SS.Dist);
DistPrime(end)       = xxx;
DistPrime            = DistPrime(:);
%--------------------------------------------------------------------------------------------------------------
% Aggregate Model
%--------------------------------------------------------------------------------------------------------------

% Capital summary
F(id.Dev.K) = log(X(id.SS.K)) - log(dot(xxx,xxx));   
% Output
F(id.Dev.Y) = log(X(id.SS.Y)) - log(X(id.SS.Z)*par.N^(1-par.alpha)*X(id.SS.K)^par.alpha); 
% TFP
F(id.Dev.Z) = log(XPrime(id.SS.Z)) - par.rho_Z*log(xxx);
% Interest rates
F(id.Dev.R) = log(X(id.SS.R)) - log(xxx);
% Consumption
F(id.Dev.C) = log(X(id.SS.C)) - log(xxx); 
% Investment
F(id.Dev.I) = log(X(id.SS.I)) - log(xxx); 

%--------------------------------------------------------------------------------------------------------------
% Individual Consumption Policies
%--------------------------------------------------------------------------------------------------------------
if par.gamma ==1
    mutil    = @(c) 1./c;  % Marginal utility
    invmutil = @(mu) 1./mu;% inverse marginal utility
else
    mutil    = @(c) 1./(c.^par.gamma); % Marginal utility
    invmutil = @(mu) 1./(mu.^(1./par.gamma)); % inverse marginal utility
end
w        = X(id.SS.Z)*par.N^(-par.alpha)*X(id.SS.K)^par.alpha*(1-par.alpha);
meshes.z = meshes.z*w;
[cpol,kprime] = EGM(xxx,X(id.SS.R),XPrime(id.SS.R),par,mpar,P,meshes,gri, mutil, invmutil);

F(id.Dev.cpol)= log(xxx) - log(xxx); 

%--------------------------------------------------------------------------------------------------------------
% Individual Consumption Policies
%--------------------------------------------------------------------------------------------------------------
Gamma = Young(kprime,gri,mpar,P);

DistNewPrime = xxx;

F(id.Dev.Dist) = DistPrime(1:end-1) - DistNewPrime(1:end-1);

end