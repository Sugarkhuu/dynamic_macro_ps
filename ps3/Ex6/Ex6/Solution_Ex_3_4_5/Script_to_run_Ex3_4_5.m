%% Solve a Bewley-Huiggett-Aiyagari model 
clear
clc
close all
addpath('Functions')
%% 1. Define parameters (Common header for Exercise 3-5)

mpar.nk   = 100;   % Number of points on the asset grid
mpar.nz   = 2;    % Number of points on the log-productivity grid
mpar.crit = 1e-5; % Numerical precision
mpar.maxk = 6;    % Maximimum assets
mpar.mink = -9/4;    % Minimum Assets (equal to Borrowing Limit)
mpar.T    = 100000;

% Economic Parameters
par.r     = 0;    % Real Rate         (for Exercise 3)
par.gamma = 2;    % Coeffcient of relative risk aversion
par.beta  = 0.95; % Discount factor
par.alpha = 0.36; % Capital Share     (for Exercise 5)
par.delta = 0.1;  % Depreciation rate (for Exercise 5)
par.b     = mpar.mink; % Borrowing Limit


%% 2. Productivity grids and transitions
prob.z  = [3/5, 2/5; 4/90,  86/90];
gri.z   = [1/9, 10/9];

%% Run the exercises (use all the same functions)
Ex3 % Bewley

Ex4 % Huggett

Ex5 % Aiyagari