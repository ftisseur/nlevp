function [coeffs,fun,F,xcoeffs] = railtrack2(n,omega)
%RAILTRACK2   Palindromic QEP from model of rail tracks.
%  [COEFFS,FUN,F,XCOEFFS] = nlevp('railtrack2',N,OMEGA) returns matrices
%  A0, A1, and A2 for the quadratic eigenvalue problem of dimension
%  DIM*705, where DIM > 1 is chosen so that DIM*705 is as close as
%  possible to N, defined by
%  Q(lam)*x = (lam^2*A2 + lam*A1 + A0)*x = 0.
%  DIM is chosen such that DIM*705 approximates N.
%  A0, A1, A2 have the form
%               |  0   0   H1  |
%          A0 = |  0   0   0   | = A2.'
%               |  0   0   0   |
%  and
%               | H0  H1.'                  |
%               | H1  H0    H1.'            |
%          A1 = |      \     \    \         | = A1.',
%               |           H1    H0  H1.'  |
%               |                 H1  H0    |
%  where H0 and H1 are 705-by-705 complex matrices depending on OMEGA.
%  The defaults are N = 51*705, OMEGA = 1000.
%  The matrices are returned in a cell array: COEFFS = {A0, A1, A2}.
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2
%  and their derivatives.
%  F is the function handle lam^2*A2 + lam*A1 + A0.
%  XCOEFFS is the cell {L0 1 L2; R0 A1 R2} that exploits the low rank of
%  A0 = L0*R0' and of A2 = L2*R2'.
%  This problem has the properties pep, qep, sparse, t-palindromic,
%  scalable, parameter-dependent, low-rank.

%  This routine is based on initial version by T.-M. Huang, July 6, 2010.

%  References
%  E. King-wah Chu, T.-M. Hwang, W.-W. Lin and C.-T. Wu, Vibration of fast trains, 
%  palindromic eigenvalue problems and structure-preserving doubling algorithms, 
%  Journal of Computational and Applied Mathematics, Vol. 219, 237-252, 2008.
%
%  T.-M. Huang, W.-W. Lin and J. Qian, Structure-preserving algorithms
%  for palindromic quadratic eigenvalue problems arising from vibration on fast 
%  trains, SIAM J. Matrix Anal. Appl., 30(4):1566–1592, 2008.
%
%  C.-H. Guo and W.-W. Lin, Solving a structured quadratic eigenvalue
%  problem by a structure-preserving doubling algorithm.
%  SIAM J. Matrix Anal. Appl., 31(5):2784–2801, 2010.

% Default values.
if nargin < 1 || isempty(n)
    n = 51*705;
else
    warning('NLEVP:truescale',['Note that the scale parameter N ',...
        'now targets the true dimension of the problem']);

end
if nargin < 2 || isempty(omega)
    omega = 1e3;
end

dim = round(n/705);
% n = dim*705;

if dim <= 1, error('N must be at least 1410.'), end

K0 = 0; K1 = 0; M0 = 0; M1 = 0; % Just to prevent lint-warnings.

load railtrack2.mat % Load K0, K1, M0, M1.

%%% Damping matrix D is constructed by r*K+(1-r)*M           %%%

r      = 2e-01;                %%% r: relaxation for damping

K0 = (K0 + K0') * 0.5;
M0 = (M0 + M0') * 0.5;
D0 = r * K0 + (1-r) * M0;
D1 = r * K1 + (1-r) * M1;
H0 = K0 + 1i * omega * D0 - omega^2 * M0;
H1 = K1' + 1i * omega * D1' - omega^2 * M1';

H0 = 0.5*(H0 + H0.');

%%% Reduced palindromic eigenvalue problem %%%

m                              = size(H0,1);

[ i_idx_H0, j_idx_H0, val_H0 ] = find( H0 );
[ i_idx_H1, j_idx_H1, val_H1 ] = find( H1 );
leng_H0                        = length(i_idx_H0);
leng_H1                        = length(i_idx_H1);
leng_A1                        = dim * leng_H0 + 2 * (dim-1) * leng_H1;

A0 = sparse(i_idx_H1, (dim-1)*m+j_idx_H1, val_H1, dim*m, dim*m );

i_idx_A1                       = zeros(leng_A1,1);
j_idx_A1                       = zeros(leng_A1,1);
val_A1                         = zeros(leng_A1,1);

i_idx_A1(1:leng_H0,1)          = i_idx_H0;
j_idx_A1(1:leng_H0,1)          = j_idx_H0;
val_A1(1:leng_H0,1)            = val_H0;

leng                           = leng_H0 + leng_H1;
i_idx_A1(leng_H0+1:leng,1)     = j_idx_H1;
j_idx_A1(leng_H0+1:leng,1)     = m + i_idx_H1;
val_A1  (leng_H0+1:leng,1)     = val_H1;

for ii = 2:dim-1
    i_idx_A1(leng+1:leng+leng_H1,1)    = (ii - 1) * m + i_idx_H1;
    j_idx_A1(leng+1:leng+leng_H1,1)    = (ii - 2) * m + j_idx_H1;
    val_A1  (leng+1:leng+leng_H1,1)    = val_H1;

    leng                               = leng + leng_H1;
    i_idx_A1(leng+1:leng+leng_H0,1)    = (ii - 1) * m + i_idx_H0;
    j_idx_A1(leng+1:leng+leng_H0,1)    = (ii - 1) * m + j_idx_H0;
    val_A1  (leng+1:leng+leng_H0,1)    = val_H0;

    leng                               = leng + leng_H0;
    i_idx_A1(leng+1:leng+leng_H1,1)    = (ii - 1) * m + j_idx_H1;
    j_idx_A1(leng+1:leng+leng_H1,1)    = ii * m + i_idx_H1;
    val_A1  (leng+1:leng+leng_H1,1)    = val_H1;

    leng                               = leng + leng_H1;
end
i_idx_A1(leng+1:leng+leng_H1,1) = (dim-1) * m + i_idx_H1;
j_idx_A1(leng+1:leng+leng_H1,1) = (dim-2) * m + j_idx_H1;
val_A1  (leng+1:leng+leng_H1,1) = val_H1;

leng                            = leng + leng_H1;
i_idx_A1(leng+1:leng+leng_H0,1) = (dim - 1) * m + i_idx_H0;
j_idx_A1(leng+1:leng+leng_H0,1) = (dim - 1) * m + j_idx_H0;
val_A1  (leng+1:leng+leng_H0,1) = val_H0;

A1 = sparse(i_idx_A1, j_idx_A1, val_A1, dim*m, dim*m );

rescale = 1d8;
A0  = A0 / rescale;
A1  = A1 / rescale;

% OUTPUT
coeffs = {A0, A1, A0.'};
fun = @(lam) nlevp_monomials(lam,2);
F = nlevp_handleQEP(coeffs);
L2 = conj(R0);
R2 = conj(L0);
xcoeffs = {L0, 1, L2; R0', A1, R2'};


end
