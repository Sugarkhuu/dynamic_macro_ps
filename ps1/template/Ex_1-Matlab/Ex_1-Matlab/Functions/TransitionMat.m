function Gamma = TransitionMat(kprime,gri,mpar,P)
% This function calculates a transition matrix from policy functions 
% and uses this to obtain the stationary distribution of income and wealth
% (Young's metheod)

[~,idk]                  = histc(kprime,gri.k); % find the next lowest point on grid for policy
idk(kprime<=gri.k(1))    = 1;           % remain in the index set
idk(kprime>=gri.k(end))  = mpar.nk-1;   % remain in the index set

% Calculate linear interpolation weights
distance                 = kprime - gri.k(idk);
weightright              = distance ./ (gri.k(idk+1)-gri.k(idk));
weightleft               = 1-weightright;

Trans_array = zeros(mpar.nk,mpar.nz,mpar.nk,mpar.nz); % Assets now x Income now x Assets next x Income next
% Fill this array with probabilities
for zz=1:mpar.nz % all current income states
    for kk=1:mpar.nk % all current asset states
        Trans_array(kk,zz,idk(kk,zz),:)   =  weightleft(kk,zz) *reshape(P(zz,:),[1 1 1 mpar.nz]); % probability to move left  to optimal policy choice
        Trans_array(kk,zz,idk(kk,zz)+1,:)   =  weightright(kk,zz)*reshape(P(zz,:),[1 1 1 mpar.nz]); % probability to move right to optimal policy choice
    end
end
Gamma  = reshape(Trans_array,[mpar.nk*mpar.nz, mpar.nk*mpar.nz]); % Turn 4-d array into a transition matrix stacking 2 dimensions
end