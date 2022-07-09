%% Exercise 6: A HANC model, Reiter's method
clear
clc
close all
addpath('Functions')
%% 1. Define parameters

%% Define Numerical Parameters
mpar.nk   = 100;  % number of points on the capital grid
mpar.nz   = 2;    % number of points on the employment grid
mpar.mink = 0;  % lowest point on the capital grid
mpar.maxk = 20;    % highest point on the capital grid
mpar.crit = 1e-10; % Precision up to which to solve the value function
mpar.T    = 10;  %length of IRF
%% Define Economic parameters
par.beta  = 0.95; % Discount factor
par.alpha = 0.36;    % Curvature of production function
par.gamma = 4;    % Coefficient of relative risk aversion
par.delta = .1;    % Depreciation
par.rho_Z = 0.75;  %persistence of technology shock
par.UIB   = 0.5;  % Unemployment benefit as fraction of employed
% unemployment budget $$U * UIB *w*h = \tau*E*w*h \implies \tau = UIB *U/E$$

par.b     = mpar.mink;

%% Produce grids
prob.z  = [3/5, 2/5; 4/90,  86/90];
aux     = prob.z^1000;
E       = aux(1,2); %Share of employed
U       = xxx; % Share of unemployed
par.N   = xxx; % aggregate labor supply, correct E for productivity h
par.tau = xxx * U/E; % tax rate to finance UIB

gri.z   = [par.UIB*10/9, 10/9 *(1-par.tau)]; % to be multiplied by w
gri.k   = exp(linspace(0,log(mpar.maxk-mpar.mink+1),mpar.nk))-1+mpar.mink; % individual capital grid

% Meshes of capital and productivity
[meshes.k,  meshes.z] = ndgrid(gri.k,gri.z);

%% Display Model
TablePar=cell2table({'Discount Factor:', par.beta; 'Returns to Scale', par.alpha; ...
    'Relative Risk Aversion', par.gamma; 'Depreciation', par.delta;...
    'Unemployment Benefit', par.UIB; 'UI-Contributions', par.tau});

disp(TablePar)
%% 2. Solve the incomplete markets model w/o aggregate risk


Kdemand         = @(R) (par.N * (par.alpha/(R+par.delta)).^(1/(1-par.alpha)));          % Calculate capital demand by firms for a given interest rate and employment
rate            = @(K) (par.alpha* par.N.^(1-par.alpha) * K.^(par.alpha-1) -par.delta); % Calculate the return on capital given K and employment N 
wage            = @(K) ((1-par.alpha)* par.N.^(-par.alpha) * K.^(par.alpha));           % Calculate the wage rate given K and employment N 
ExcessDemand    = @(K) (K_Agg(rate(K),wage(K),par,mpar,prob.z,meshes,gri) - K);     % Calculate the difference between capital supply and demand for wages and returns given by assumed capital demand

Rstar           = rate(fzero(ExcessDemand,[Kdemand(0.052),Kdemand(-0.01)])); % find equilibrium amount of capital (and corresponding rate)
Kstar           = Kdemand(Rstar);
[~,saving_star,marginal_k, StDist, Gamma,c_star] = K_Agg(rate(Kstar),wage(Kstar),par,mpar,prob.z,meshes,gri);
%% 3. introduce aggregate Risk: Reiter method
%Create indexes to access entries in steady state and deviation threof
%vector
count        = 0;
% States
id.SS.Dist   = count+(1:(numel(StDist)-1));    count   = count + length(id.SS.Dist); % leave out one element from distribution
id.SS.Z      = count+1;                        count   = count + 1;
% Controls
id.SS.cpol   = count+(1:numel(c_star));        count   = count + length(id.SS.cpol);
id.SS.R      = count+1;                        count   = count + 1;
id.SS.K      = count+1;                        count   = count + 1;
id.SS.Y      = count+1;                        count   = count + 1;
id.SS.C      = count+1;                        count   = count + 1;
id.SS.I      = count+1;                        count   = count + 1;
id.Dev       = id.SS; % for Reiter no difference between deviation and steady state size (no compression)

XSS             = zeros(count,1);
XSS(id.SS.cpol) = xxx;
XSS(id.SS.R)    = 1+Rstar;
XSS(id.SS.K)    = Kstar;
XSS(id.SS.Y)    = par.N^(1-par.alpha)*Kstar^par.alpha;
XSS(id.SS.C)    = par.N^(1-par.alpha)*Kstar^par.alpha - par.delta*Kstar;
XSS(id.SS.I)    = par.delta*Kstar;
XSS(id.SS.Z)    = 1;
XSS             = log(XSS);
XSS(id.SS.Dist) = StDist(1:end-1);

mpar.nstates    = length(id.Dev.Dist) + 1;
mpar.ncontrols  = length(XSS)-mpar.nstates;
[A,B]           = Derivatives(XSS,  mpar, par, id, gri, meshes,prob.z);

[S2C,S2SPrime]  = Klein(A,B,mpar);


%% 4. Plot IRF
s_t = zeros(mpar.nstates,mpar.T+1);
x_t = zeros(mpar.nstates+mpar.ncontrols,mpar.T);
s_t(id.Dev.Z,1) = 1;
for t=1:mpar.T
    s_t(:,t+1) = xxx*s_t(:,t);
    x_t(:,t) = [eye(mpar.nstates) ; xxx] *s_t(:,t) ;
end

plot(x_t([id.Dev.Y, id.Dev.C, id.Dev.I, id.Dev.R],:)')
legend({'Output', 'Consumption', 'Investment', 'interest rate'},'Location', 'northeast')
