function [coeffs,fun,F,xcoeffs] = shaft()
%SHAFT  QEP from model of a shaft on bearing supports with a damper.
%  [COEFFS,FUN,F,XCOEFFS] = nlevp('shaft') generates the coefficient matrices of a
%  quadratic matrix polynomial lambda^2*A2 + lambda*A1 + A0 arising from
%  a finite element model of a shaft on bearing supports with a damper.
%  This problem is real symmetric and has dimension 400.
%  The matrices are sparse and returned in a cell array:
%  COEFFS = {A0, A1, A2}.
%  F is the function handle lambda^2*A2 + lambda*A1 + A0.
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2
%  and their derivatives.
%  XCOEFFS is the cell {1, l1, 1; A0, r1', A2} to exploit the low rank of
%  A1 = l1*r1'.
%  This problem has the properties pep, qep, real, symmetric, sparse,
%  banded (6 bands), low-rank.

load('shaft.mat');


coeffs = {K, C, M};
fun = @(lam)nlevp_monomials(lam,2);
F = nlevp_handleQEP(coeffs);
% Building xcoeffs
en = sparse(zeros(400,1));
en(20,1) = 1;
a = C(20,20);
xcoeffs = {1, a*en, 1; K, en', M};
end
