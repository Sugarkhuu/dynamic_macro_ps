
%% I. Exercise 3: Bewley model
% Specific parameter adjustments
mpar.mink = 0;
mpar.maxk = 3;    % Maximimum assets
par.b     = mpar.mink;
disp('-------------------------')
disp('Exercise 3: Bewley model')
disp('-------------------------')
disp('Economic parameters')
disp(par) % Display economic parameters
%% 1. Generate grids, Meshes and Income
gri.k   = exp(linspace(0,log(mpar.maxk-mpar.mink+1),mpar.nk))-1+mpar.mink; %Define asset grid on log-linearspaced

% Meshes of capital and productivity
[meshes.k,  meshes.z] = ndgrid(gri.k,gri.z);

%% 2. Equilibrium in Bewley model with young method
[K,kprime,marginal_k, ~, ~] = K_Agg(0,1,par,mpar,prob.z,meshes,gri); % Solve for equlibrium

%% 3. Plot Policy Functions
figure(1)
plot(gri.k,kprime) % Plot policy functions
hold on
plot(gri.k,gri.k,'k--') % Add 45° line
hold off
title('Policy Value Functions') % Title and legend of the graph
legend({'low productivity','high productivity', '45 degree'},'Location','northwest')
xlabel('assets')
ylabel('consumption')
%% 4. Simulate the economy
% Define an off-grid savings function by interpolation
tic
Saving  = griddedInterpolant({gri.k,gri.z},kprime); % Continuous Savings Function as Linear Interpolant
% Simulate the exogeneous state
PI      = cumsum(prob.z,2);         % Cummulative Transition matrix
epsilon = rand(1,mpar.T);           % Random numbers for simulation
S       = randi(mpar.nz,1,mpar.T);  % Starting value
k       = zeros(1,mpar.T);          % Starting value for NEXT PERIODS ASSETS assets
for t=2:mpar.T
    S(t)    = sum(PI(S(t-1),:)<epsilon(t))+1;   % Update productivity state
    k(t)    = Saving({k(t-1),gri.z(S(t))});     % Update Assest Holdings
end
time(1)     = toc;

% Plot histogram
figure(2)
histogram(k(10001:end),'BinEdges',(gri.k(1:end-1) + gri.k(2:end))/2,'Normalization','probability')
title('Distribution of asset holdings')
hold on

%% 5. Compare solutions
tic % time young method 
[Gamma, StDist] = Young(kprime,gri,mpar,prob.z);
time(2) = toc;

figure(2)
bar(gri.k,marginal_k)
title('Distribution of asset holdings')
legend('Simulation','Direct Calculation')
hold off

disp('Average Asset holdings')
disp('----------------------')
disp(['Simulation: ' num2str(mean(k(10001:end)))])
disp(['Young     : ' num2str(dot(marginal_k',gri.k))]);
disp('----------------------')
disp('Time for calculation')
disp('----------------------')
disp(['Simulation: ' num2str(time(1))])
disp(['Young     : ' num2str(time(2))]);
