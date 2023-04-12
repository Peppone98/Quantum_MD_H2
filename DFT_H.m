%% SCRIPT FOR H atom DFT

g = 100; 
g3 = g^3;
p = linspace(-5, 5, g);
[X, Y, Z] = meshgrid(p, p, p);
h = p(2) - p(1);
X = X(:);
Y = Y(:);
Z = Z(:);
R = sqrt(X.^2 + Y.^2 + Z.^2);
%R_prime = sqrt((X).^2 + Y.^2 + Z.^2);
%Vext = -1./R - 1./R_prime;
Vext = -1./R;
e = ones(g, 1);
L = spdiags([e -2*e e], -1:1, g, g)/h^2;
I = speye(g);
L3 = kron(kron(L, I), I) + kron(kron(I, L), I) + kron(kron(I, I), L);
E = eigs(-0.5*L3 + spdiags(Vext, 0, g3, g3), 1, 'sa');
disp(['Total energy for H: ' num2str(E) ' H']);

 
 