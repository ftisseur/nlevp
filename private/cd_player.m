function [coeffs,fun,F] = cd_player
%CD_PLAYER      QEP from model of CD player.
%  [COEFFS,FUN,F] = nlevp('cd_player') constructs a 60-by-60 quadratic matrix
%  polynomial lambda^2*M + lambda*D + K arising in the
%  study of a CD player control task.
%  The matrices are returned in a cell array: COEFFS = {K, D, M}.
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2
%  and their derivatives.
%  F is the function handle K + lambda*D + lambda^2*M.
%  This problem has the properties pep, qep, real.

%  Reference:
%  Y. Chahlaoui and P. M. Van Dooren. A collection of benchmark examples
%  for model reduction of linear time invariant dynamical systems.
%  MIMS EPrint 2008.22, Manchester Institute for Mathematical Sciences,
%  The University of Manchester, UK, 2002.

load cd_player.mat

n = size(D,1);
M = eye(n);

coeffs = {K, D, M};
fun = @(lam) nlevp_monomials(lam,2);
F =  nlevp_handleQEP(coeffs);

end
