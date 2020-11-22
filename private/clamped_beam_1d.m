function [coeffs,fun,F,xcoeffs] = clamped_beam_1d(n,tau)
%CLAMPED_BEAM_1D  NEP from 1D clamped beam model with delayed feedback control.
%  [COEFFS,FUN,F] = nlevp('clamped_beam_1D',N,tau) returns N-by-N
%  matrices I, A0 and A1 for the nonlinear eigenvalue problem defined by
%  [lambda*I+A0+exp(-tau*lambda)*A1]*x = 0, where I is the identity matrix,
%  A0 is tridiagonal and A1 = -e_n^T*e_n, where e_n is the nth column of
%  the identity matrix.
%  The matrices are returned in a cell array: COEFFS = {EYE(N),A0,A1}.
%  FUN is a function handle to evaluate the functions lam,1, and exp(lambda)
%  and their first derivatives.
%  F is the function handle lambda*I+A0+exp(-tau*lambda)*A1.
%  XCOEFFS return the cell {1, 1, -en; EYE(N), A0, en'} to exploit the
%  low-rank of A0.
%  This problem has the properties nep, real, sparse, parameter_dependent,
%  scalable, tridiagonal, banded, low-rank.

%   Reference:   Sec. 5.2 in
%   R. Van Beeumen, Elias Jarlebring and W. Michiels.
%   A rank-exploiting infinite Arnoldi algorithm for nonlinear
%   eigenvalue problems. Numer. Linear Algebra Appl., 23:607-628, 2016.


if nargin < 1 || isempty(n), n = 10001; end;
if nargin < 2 || isempty(tau), tau = 1; end;

A0 = spdiags(ones(n,1)*[-1,2,-1],-1:1,n,n); A0(n,n) = -n; A0(n,n-1) = n;
A1 = sparse(n,n,-1,n,n);

coeffs = {speye(n), A0, A1};
fun = @(lam) clamped_beam_1D_fun(tau, lam);
F = @(lam) lam*coeffs{1} + coeffs{2} + coeffs{3}*exp(-tau*lam);

en = sparse(n,1); 
en(n) = 1;
xcoeffs1 = {1, 1, -en};
xcoeffs2 = {coeffs{1}, coeffs{2}, en'};
xcoeffs = {xcoeffs1{:}; xcoeffs2{:}};
end

function varargout = clamped_beam_1D_fun(tau, lam)
lam = lam(:);
n = length(lam);
varargout{1} = [lam,ones(n,1),exp(-tau*lam)];

if nargout >= 2
    f = -tau*exp(-tau*lam);
    varargout{2} = [ones(n,1), zeros(n,1), f];
    for i=2:nargout-1
        f = -tau*f;
        varargout{i+1} = [zeros(n,2),f];
    end
end
end


