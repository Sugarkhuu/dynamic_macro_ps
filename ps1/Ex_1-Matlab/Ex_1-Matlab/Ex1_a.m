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
% attention logspace gives logspaced grid on a basis 10
% use exp(linspace(...)) instead
% make the adjustment to move the log-spaced grid to cover negative asset
% holdings, too.
gri.k   = exp(linspace(0,log(),mpar.nk)) - xxx; %Define asset grid on log-linearspaced
prob.z  = [3/5, 2/5; 4/90,  86/90];
gri.z   = [1/9, 10/9];


%% 3. Define utility functions

if par.gamma ==1
    util     = @(c) xxx; % Utility
    mutil    = @(c) xxx;  % Marginal utility
    invmutil = @(mu) xxx;% inverse marginal utility
else
    util     = @(c) xxx; % Utility
    mutil    = @(c) xxx; % Marginal utility
    invmutil = @(mu) xxx; % inverse marginal utility
end


%% 4b. Value Function Iteration (on-grid)
% Meshes of capital and productivity
[meshes.kprime,   meshes.k, meshes.z] = ndgrid(xxx,xxx,xxx);
Y               = xxx; % Cash at hand (Labor income plus assets cum dividend)
tic % Start timer
V               = zeros(mpar.nk,mpar.nz); % Initialize Value Function
distVF_on       = 1; % Initialize Distance
iterVF          = 1; % Initialize Iteration count
Chat            = xxx; % Consumption: Cash at hand minus investment
U               = util(Chat); % evaluate utility
U(Chat<=0)      = - 1.0e10;  % replace negative consumption utility by high negative value
while distVF_on(iterVF)>mpar.crit % Value Function iteration loop: until distance is smaller than crit.
    % Update Value Function using on-grid search
    EV              = xxx;           
    Vhat            = xxx + repmat(reshape(EV,[mpar.nk,1,mpar.nz]), [1,mpar.nk,1]); 
    [Vaux,pol_ind]  = max(xxx,[],1); % Optimize given cont' value
    Vnew            = squeeze(Vaux); 
    kprime_on       = squeeze(gri.k(xxx));
    dd              = max(abs(Vnew(:)-V(:))); % Calculate distance between old guess and update
    V               = Vnew; % Update Value Function
    iterVF          = iterVF+1; %Count iterations
    distVF_on(iterVF)  = dd;   % Save distance
end
V_on       = reshape(V,[mpar.nk,mpar.nz]);
time(1)    = toc; % Save Time used for VFI


%% 4b. Value Function Iteration (off-grid)
% Meshes of capital and productivity
[meshes.k,  meshes.z] = ndgrid(gri.k,gri.z);
Y = meshes.z + meshes.k*(1+par.r); % Cash at hand (Labor income plus assets cum dividend)
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
V      = reshape(V,[mpar.nk,mpar.nz]);
time(2)=toc; % Save Time used for VFI
%% 5. Plot Value and Policy Functions from VFI

figure(1)
plot(gri.k,V_on,'--') % Plot policy functions
hold on
plot(gri.k,V) % Plot policy functions
hold off
title('Value Function from VFI') % Title and legend of the graph
legend({'low productivity (on)','high productivity (on)','low productivity (off)','high productivity (off)'},'Location','northwest')
xlabel('assets')



figure(2)
plot(gri.k,kprime_on,'--') % Plot policy functions
hold on
plot(gri.k,kprime) % Plot policy functions
plot(gri.k,gri.k,'k--') % Add 45° line
hold off
title('Policy Function from VFI') % Title and legend of the graph
legend({'low productivity (on)','high productivity (on)','low productivity (off)','high productivity (off)','45 degree'},'Location','northwest')
xlabel('assets')
ylabel('saving')
%% 6. Display time and convergence stats
figure(3)
disp('Time for solution')
disp(['off grid: ', num2str(time(1))])
disp(['on  grid: ', num2str(time(2))])

semilogy(distVF_on(2:end))
hold on
semilogy(distVF(2:end))
hold off
title('distance after iterations')
legend({'on grid', 'off grid'})




