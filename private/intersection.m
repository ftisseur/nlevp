function [coeffs,fun,F] = intersection
%INTERSECTION 10-by-10 QEP from intersection of three surfaces.
%  [COEFFS,FUN,F] = nlevp('intersection') generates a 10x10 quadratic
%  matrix polynomial Q(lambda) = lambda^2*A2 + lambda*A1 + A0 that arises
%  in the problem of finding the intersection between a cylinder, a sphere
%  and a plane described by the equations
%      1.6e-3*x^2 + 1.6e-3*y^2 = 1,
%      5.3e-4*x^2 + 5.3e-4*y^2 + 5.3e-4*z^2 + 2.7e-2*x = 1,
%      -1.4e-4*x + 1.0e-4*y + z = 3.4e-3.
%  The quadratic Q(lambda) is obtained from the Macaulay resultant and
%  has two real eigenpairs (lambda,v). With the normalization u(10) = 1,
%  (x,y,z) = (lambda, u(8), u(9)) is a solution to the above equations.
%  The matrices are returned in a cell array: COEFFS = {A0, A1, A2}.
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2.
%  F is the function handle Q(lambda) = lambda^2*A2 + lambda*A1 + A0.
%  This problem has the properties pep, qep, real.

%  Reference:
%  D. Manocha, Solving Systems of Polynomial Equations, IEEE Computer
%  Graphics and Applications, 14(2):46-55, 1994.

i0 = [1 4 2 5 7 3 6 8 4 7 5 6 9 7 8 9 1 4 8 10 2 5 9 10 3 6 10];
j0 = [1 1 2 2 2 3 3 3 4 4 5 6 6 7 7 7 8 8 8 8 9 9 9 9 10 10 10];
a0 = [1.6e-3 5.3e-4 1.6e-3 5.3e-4 1.0e-4 1.6e-3 5.3e-4 1.0e-4 5.3e-4 1 ...
      5.3e-4 5.3e-4 1 -3.4e-3 1 1.0e-4 -1 -1 -3.4e-3 1.0e-4 -1 -1  ...
      -3.4e-3 1 -1 -1 -3.4e-3];
A0 = full(sparse(i0,j0,a0));


i1 = [7 4 8 5 9  6 10];
j1 = [7 8 8 9 9 10 10];
a1 = [-1.4e-4 0.027 -1.4e-4 0.027 -1.4e-4 0.027 -1.4e-4];
A1 = full(sparse(i1,j1,a1));


i2 = [1 4 2 5 3   6 10];
j2 = [8 8 9 9 10 10 10];
a2 = [1.6e-3 5.3e-4 1.6e-3 5.3e-4 1.6e-3 5.3e-4 0];
A2 = full(sparse(i2,j2,a2));

coeffs = {A0, A1, A2};
fun = @(lam)nlevp_monomials(lam,2);
F = nlevp_handleQEP(coeffs);

end
