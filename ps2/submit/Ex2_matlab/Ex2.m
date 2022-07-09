%% Solve the Consumption Savings Model using an Endogenous Grid Method (EGM)

clear
clc
close all
addpath('Functions')
%% 1. Define parameters

% Numerical parameters
mpar.nk   = 100;   % Number of points on the asset grid
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
% Meshes of capital and productivity
[meshes.k,  meshes.z] = ndgrid(gri.k,gri.z);
Y = meshes.z + meshes.k*(1+par.r); % Cash at hand (Labor income plus assets cum dividend)

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

%% 4. Value Function Iteration

tic % Start timer
V    = zeros(mpar.nk,mpar.nz); % Initialize Value Function
distVF = 1; % Initialize Distance
iterVF = 1; % Initialize Iteration count
while distVF(iterVF)>mpar.crit % Value Function iteration loop: until distance is smaller than crit.
    % Update Value Function using off-grid search
    [Vnew,kprime] = VFI_update_spline(V,Y,util,par,mpar,gri,prob); % Optimize given cont' value
    dd            = max(abs(Vnew(:)-V(:))); % Calculate distance between old guess and update

    V             = Vnew; % Update Value Function
    iterVF        = iterVF+1; %Count iterations
    distVF(iterVF)= dd;   % Save distance
end
time(1)=toc; % Save Time used for VFI
%% 5. Plot Policy Functions from VFI

figure(1)
plot(gri.k,kprime) % Plot policy functions
hold on
plot(gri.k,gri.k,'k--') % Add 45° line
hold off
title('Policy Function from VFI') % Title and legend of the graph
legend({'low productivity','high productivity', '45 degree'},'Location','northwest')
xlabel('assets')
ylabel('saving')

%% 6. Endogenous Grid method using linear interpolation
% $$\frac{\partial u}{\partial c}\left[C^*(k',z)\right]=(1+r) \beta E_{z}\left\{\frac{\partialu}{\partial
% c}\left[C(k',z')\right]\right\}$$

tic % Reset timer
C       = meshes.z + meshes.k*(1+par.r); %Initial guess for consumption policy: roll over assets
C(C<0.01) = 0.01; % to avoid negative consumption
Cold    = C; % Save old policy
distEG  = 1; % Initialize Distance
iterEG  = 1; % Initialize Iteration count
while distEG(iterEG)>mpar.crit
    C      = EGM(Cold, mutil,invmutil,par,mpar,prob.z,meshes,gri); % Update consumption policy by EGM
    dd     = max(abs(C(:) - Cold(:))); % Calculate Distance

    Cold   = C; % Replace old policy
    iterEG = iterEG+1; %count iterations
    distEG(iterEG) = dd;
end
[C,Kprimestar] = EGM(C,mutil,invmutil,par,mpar,prob.z,meshes,gri);
time(2)        = toc; %Time to solve using EGM
%% 7. Plot Policy Functions from Collocation and compare to VFI

figure(2) %Plot Policy Functions from Collocation
plot(gri.k,Kprimestar) % Plot Policy function from Collocation
hold on
plot(gri.k,gri.k,'k--') % Add 45° line
hold off
title('Policy Function from EGM') % Title and Legend
legend({'low productivity','high productivity', '45 degree'},'Location','northwest')
xlabel('assets')
ylabel('saving')

figure(3) %Plot Differences in Policies
plot(gri.k,Kprimestar - kprime)
title('Difference in Policy Function')

%% 8. Compare times of algorithms
disp('Time to solve (VFI, EGM)')
disp(time)
disp('Iterations to solve (VFI, EGM)')
disp([iterVF iterEG])



