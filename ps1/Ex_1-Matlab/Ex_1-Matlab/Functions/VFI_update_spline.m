%% VF Update
function [Vnew,kprime] = VFI_update_spline(V,Y,util,par,mpar,gri,prob)
  % VFI_update_spline updates the value function (one VFI iteration) for the
  %  consumption-savings problem.
  % V (dimensions: k x z) is the old value function guess.
  % Y (dimensions: k x z) is a matrix of cash at hand UTIL is the felicity function.
  % PAR and MPAR are structures containing economic and numerical  parameters.
  % PROB (dimensions: z x z') is the transition probability matrix.

V      = reshape(V,[mpar.nk,mpar.nz]); % make sure that V has the right format dim1: k, dim2:z
kprime = zeros(size(V)); % allocate policy matrix
Vnew   = zeros(size(V)); % allocate new value matrix
EV     = xxx;   % Calculate expected continuation value

for zz=1:mpar.nz % loop over Incomes
    ev_int = griddedInterpolant({xxx},xxx,'spline'); % Prepare interpolant
    for kk=1:mpar.nk % loop of Assets
        f             = @(k)(-util(xxx)-ev_int({k})); % Define function to be minimized
        [kp,v]        = fminbnd(f,xxx,xxx); % Find minimum of f for savings between borrowing constraint and cash at hand 
        Vnew(kk,zz)   = xxx;  % Save value
        kprime(kk,zz) = xxx; % Save policy
    end
end
Vnew=Vnew(:);
end