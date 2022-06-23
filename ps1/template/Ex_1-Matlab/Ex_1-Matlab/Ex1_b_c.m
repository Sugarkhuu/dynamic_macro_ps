%% Solve the Consumption Savings Model using Value Function Iteration

clear
clc
close all
addpath('Functions')
%% 1. Define parameters

% Numerical parameters
mpar.nk   = 30;   % Number of points on the asset grid
mpar.nz   = 2;    % Number of points on the log-productivity grid
mpar.crit = 1e-5; % Numerical precision
mpar.maxk = 6;    % Maximimum assets
mpar.mink = -9/4;    % Minimum Assets (equal to Borrowing Limit)
disp('Numerical parameters')
disp(mpar) % Display numerical parameters
% Economic Parameters
par.r     = 4/90;% Real Rate
par.gamma = 1;    % Coeffcient of relative risk aversion
par.beta  = 0.95; % Discount factor
par.b     = mpar.mink; % Borrowing Limit
disp('Economic parameters')
disp(par) % Display economic parameters

%% 2. Generate grids, Meshes and Income
gri.k   = exp(linspace(log(1),log(mpar.maxk-mpar.mink+1),mpar.nk))-1+mpar.mink; %Define asset grid on log-linearspaced
prob.z  = [3/5, 2/5; 4/90,  86/90];
gri.z   = [1/9, 10/9];


%% 3. Define utility functions

if par.gamma ==1
    util     = @(c)log(c); % Utility
    mutil    = @(c) 1./c;  % Marginal utility
    invmutil = @(mu) 1./mu;% inverse marginal utility
else
    util     = @(c) 1/(1-par.gamma).*c.^(1-par.gamma); % Utility
    mutil    = @(c) 1./(c.^par.gamma); % Marginal utility
    invmutil = @(mu) 1./(mu.^(1./par.gamma)); % inverse marginal utility
end


%% 4a. Value Function Iteration (off-grid)
% Meshes of capital and productivity
[meshes.k,  meshes.z] = ndgrid(gri.k,gri.z);
Y = meshes.z + meshes.k*(1+par.r); % Cash at hand (Labor income plus assets cum dividend)
tic % Start timer
V    = zeros(mpar.nk,mpar.nz); % Initialize Value Function
distVF = 1; % Initialize Distance
iterVF = 1; % Initialize Iteration count
while distVF(iterVF)>mpar.crit % Value Function iteration loop: until distance is smaller than crit.
    % Update Value Function using off-grid search
    [Vnew,~]      = VFI_update_spline(V,Y,util,par,mpar,gri,prob); % Optimize given cont' value
    dd            = max(abs(Vnew(:)-V(:))); % Calculate distance between old guess and update

    V             = Vnew; % Update Value Function
    iterVF        = iterVF+1; %Count iterations
    distVF(iterVF)= dd;   % Save distance
end
V      = reshape(V,[mpar.nk,mpar.nz]);
time(1)=toc; % Save Time used for VFI
keyboard
%% 4b. Value Function Iteration (Broyden)
% Meshes of capital and productivity
tic % Start timer

Dist_V         = @(V)  (V(:)- reshape(VFI_update_spline(V,Y,util,par,mpar,gri,prob),[mpar.nk*mpar.nz,1]));
V              = zeros(mpar.nk,mpar.nz); % Initialize Value Function
[V,~,~,dist_B] = broyden(Dist_V,V(:),mpar.crit,1e-14,250);
[~,kprime]     = VFI_update_spline(V,Y,util,par,mpar,gri,prob); % Optimize given cont' value
V              = reshape(V,[mpar.nk,mpar.nz]);
time(2)=toc; % Save Time used for VFI

%% 4a. Policy Function Iteration (off-grid)
% Meshes of capital and productivity
tic % Start timer
V    = zeros(mpar.nk,mpar.nz); % Initialize Value Function
dist_PFI = 1; % Initialize Distance
iterPF = 1; % Initialize Iteration count
while dist_PFI(iterPF)>mpar.crit % Value Function iteration loop: until distance is smaller than crit.
    % Update Value Function using off-grid search
    [~,kprime] = VFI_update_spline(V,Y,util,par,mpar,gri,prob); % Optimize given cont' value
    % Value Function update
    U          = util(Y-kprime); % Payoff under optimal policy
    Gamma      = TransitionMat(kprime,gri,mpar,prob.z); % Transition matrix
    Vnew       = (eye(numel(V(:))) - par.beta*Gamma)\U(:);
    dd         = max(abs(Vnew(:)-V(:))); % Calculate distance between old guess and update
    V          = Vnew; % Update Value Function
    iterPF     = iterPF+1; %Count iterations
    dist_PFI(iterPF)= dd;   % Save distance
end
V      = reshape(V,[mpar.nk,mpar.nz]);
time(3)=toc; % Save Time used for VFI
[~,kprime_PFI]     = VFI_update_spline(V,Y,util,par,mpar,gri,prob); 
%%Plot policies
figure(1)
plot(gri.k,kprime) % Plot policy functions
hold on
plot(gri.k,kprime_PFI,'--') % Plot policy functions
plot(gri.k,gri.k,'k--') % Add 45° line
hold off
title('Policy Function from VFI') % Title and legend of the graph
legend({'low productivity (VFI)','high productivity (VFI)','low productivity (PFI)','high productivity (PFI)','45 degree'},'Location','northwest')
xlabel('assets')
ylabel('saving')
%% 6. Display time and convergence stats
figure(2)
disp('Time for solution')
disp(['off grid VFI: ', num2str(time(1))])
disp(['Broyden off  grid: ', num2str(time(2))])
disp(['PFI off  grid: ', num2str(time(3))])

semilogy(distVF(2:end))
hold on
semilogy(dist_B(1:end))
semilogy(dist_PFI(2:end))
hold off
title('distance after iterations')
legend({'VFI', 'Broyden','PFI'})




