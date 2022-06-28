
%% II.) Exercise 4: Huggett model
% Specific parameter adjustments
mpar.maxk = 10;         % Maximimum assets
mpar.mink = -9/4;       % Minimum Assets (equal to Borrowing Limit)
par.b     = mpar.mink;  % Borrowing Limit
par       = rmfield(par,'r');    
disp('-------------------------')
disp('Exercise 4: Hugget model')
disp('-------------------------')
disp('Economic parameters')
disp(par) % Display economic parameters

%% 1. Generate grids, Meshes and Income
gri.k   = exp(linspace(0,log(mpar.maxk-mpar.mink+1),mpar.nk))-1+mpar.mink; % Define asset grid on log-linearspaced
% Meshes of capital and productivity
[meshes.k,  meshes.z] = ndgrid(gri.k,gri.z);
%% 2. Calculate Excess demand
ExcessDemand  = @(R) (K_Agg(R,1,par,mpar,prob.z,meshes,gri));   % Calculate the difference between capital supply and demand for wages and returns given by assumed capital demand

Rgrid = -0.01:.001:.049;                                        % a grid for interest rates for plotting
ExD   = arrayfun(ExcessDemand,Rgrid);                           % calculate excess demand for these rates
%% 3. Find equilibrium
Rstar_Huggett = fzero(ExcessDemand,[-0.0,0.049]);               % Calculate equilibrium rate
%% 4. Plot
figure(3)
plot(ExD,Rgrid,'LineWidth',2)
hold on
plot(0*Rgrid,Rgrid,'k')
plot([mpar.mink, mpar.maxk]*0.8,[Rstar_Huggett, Rstar_Huggett], 'k:')
hold off
xlabel('Funds')
ylabel('interest rate')
legend({'Supply of Funds', 'Demand for funds','Equilibrium Rate'},'Location','southeast')
disp('Equilibrium Interest Rate')
disp(Rstar_Huggett)



