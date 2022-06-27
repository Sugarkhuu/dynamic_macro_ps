%% III.) Exercise 5: Aiyagari Model
% Specific parameter adjustments
mpar.maxk = 20;    % Maximimum assets
mpar.mink = -9/4;    % Minimum Assets (equal to Borrowing Limit)
par.b     = mpar.mink; % Borrowing Limit
par.gamma = 4;

disp('--------------------------')
disp('Exercise 5: Aiyagari model')
disp('--------------------------')
disp('Economic parameters')
disp(par) % Display economic parameters

%% 1. Generate grids, Meshes and Income
gri.k   = exp(linspace(0,log(mpar.maxk-mpar.mink+1),mpar.nk))-1+mpar.mink; %Define asset grid on log-linearspaced
% Meshes of capital and productivity
[meshes.k,  meshes.z] = ndgrid(gri.k,gri.z);
% Calculate stationary labor supply
aux = prob.z^1000;
N   = dot(xxx,xxx);

%% 2. Calculate Excess demand
Kdemand         = @(R) (xxx);          % Calculate capital demand by firms for a given interest rate and employment
rate            = @(K) (xxx); % Calculate the return on capital given K and employment N 
wage            = @(K) (xxx);           % Calculate the wage rate given K and employment N 
ExcessDemand    = @(K) (xxx);     % Calculate the difference between capital supply and demand for wages and returns given by assumed capital demand


Rgrid = -0.01:.0025:.045;           % a grid for interest rates for plotting
KD    = arrayfun(Kdemand, Rgrid);   % calculate capital demand for these rates
ExD   = arrayfun(ExcessDemand,KD);  % calculate excess demand for these amounts of capital

%% 3. Find equilibrium
Rstar_Aiyagari  = rate(fzero(ExcessDemand,[Kdemand(0.045),Kdemand(0.00)])); % find equilibrium amount of capital (and corresponding rate)
disp('Equilibrium Interest Rate')
disp(Rstar_Aiyagari)

%% 4. Plot
figure(4)
plot(ExD+KD,Rgrid,'LineWidth',2)
hold on
plot(KD,Rgrid,'k')
plot([mpar.mink, mpar.maxk]*0.8,[Rstar_Aiyagari, Rstar_Aiyagari], 'k:')
hold off
xlabel('Funds')
ylabel('interest rate')
legend({'Supply of Funds', 'Demand for funds','Equilibrium Rate'},'Location','southeast')





