function [A,B] = Derivatives(XSS,  mpar, par, id, gri, meshes,P)
nx = length(XSS);
F  = @(XhatPrime,Xhat) Fsys_Reiter(Xhat, XhatPrime, XSS, mpar, par, id,gri,meshes,P);
F0 = F(zeros(nx,1),zeros(nx,1));

epsilon = 0.0001;

A = zeros(nx);
for j=1:nx
    h = zeros(nx,1);
    h(j) = epsilon;
    A(:,j) = (F(h , xxx ) - F0)/epsilon;
end

B = zeros(nx);
for j=1:nx
    h = zeros(nx,1);
    h(j) = epsilon;
    B(:,j) = (F(xxx,h ) - F0)/epsilon;
end

end