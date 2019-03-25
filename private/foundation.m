function [coeffs,fun,F] = foundation
%FOUNDATION  Sparse QEP from model of machine foundations.
%  [COEFFS,FUN,F] = nlevp('foundation') generates the coefficient matrices of
%  a quadratic matrix polynomial lambda^2*C + lambda*B + A arising in a
%  model of reinforced concrete machine foundations resting on the ground.
%  The model includes hysteretic damping.
%  This problem is complex symmetric and has dimension 3627.
%  The matrices are sparse and returned in a
%  cell array: COEFFS = {A, B, C}.
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2
%  and their derivatives.
%  F is the function handle A + lambda*B + lambda^2*C.
%  This problem has the properties pep, qep, symmetric, sparse.

%  Reference:
%  A. Feriani, F. Perotti and V. Simoncini,
%  Iterative system solvers for the frequency analysis of
%  linear mechanical systems,
%  Computer Methods in Appl. Mech. Eng., 190 (2000), pp. 1719-1739.

load('foundation.mat')

coeffs = {K,D,M};
fun = @(lam) nlevp_monomials(lam,2);
F = nlevp_handleQEP(coeffs);

end
