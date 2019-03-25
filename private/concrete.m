function [coeffs,fun,F,xcoeffs] = concrete(mu)
%CONCRETE  Sparse QEP from model of a concrete structure.
%  [COEFFS,FUN,F,XCOEFFS] = nlevp('concrete',MU) generates the coefficient
%  matrices of a quadratic matrix polynomial lambda^2*C + lambda*B + A
%  arising in a model of a concrete structure supporting a machine
%  assembly. 
%  This problem is complex symmetric and has dimension 2472.
%  The matrices are sparse and returned in a
%  cell array: COEFFS = {A, B, C}.
%  C = mass matrix, real diagonal, low rank.
%  B = viscous damping matrix, purely imaginary
%      and diagonal, B = i*C_v = L*R', low rank.
%  A = stiffness + uniform hysteretic damping
%      matrix, A = (1 + i*MU)*K.
%  By default, MU = 0.04.
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2
%  and their derivatives.
%  F is the function handle A + lambda*B + lambda^2*C.
%  XCOEFFS returns the cell {1, L, 1; A, R', C} to exploit the low rank of
%  the problem.
%  This problem has the properties pep, qep, symmetric,
%  parameter-dependent, sparse, low-rank.

%  Reference:
%  A. Feriani, F. Perotti and V. Simoncini,
%  Iterative system solvers for the frequency analysis of
%  linear mechanical systems,
%  Computer Methods in Appl. Mech. Eng., 190 (2000), pp. 1719-1739.

if nargin < 1
    mu = 0.04; 
end

load('concrete.mat');  % loads K_c1, C_c1, M_c1

coeffs = {(1+mu*1i)*K_c1, C_c1, M_c1};
fun = @(lam) nlevp_monomials(lam,2);
F =  nlevp_handleQEP(coeffs);
xcoeffs = {1 L 1; coeffs{1} R coeffs{3}};
end
