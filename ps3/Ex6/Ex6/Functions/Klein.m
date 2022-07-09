function [S2C,S2Sprime]=Klein(A,B,mpar)
[s, t, Q, Z] = qz(A,-B);

relev = abs(diag(s))./abs(diag(t));
ll    = sort(relev);
slt   = relev>=1;
nk    = sum(slt); 

if nk>mpar.nstates
    warning(['The Equilibrium is Locally Indeterminate!' ])
elseif nk<mpar.nstates
    warning(['No Local Equilibrium Exists!'])
end
[s,t,~,Z] = ordqz(s,t,Q,Z,slt);

z21=Z(nk+1:end,1:nk);
z11=Z(1:nk,1:nk);
s11=s(1:nk,1:nk);
t11=t(1:nk,1:nk);

%Checks
if rank(z11)<nk
    warning('invertibility condition violated')
end
z11i=z11\eye(nk);
S2C=real(z21*z11i); %States2Controls
S2Sprime=real(z11*(s11\t11)*z11i); %LOM states
end