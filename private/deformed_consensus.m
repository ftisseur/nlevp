function [coeffs,fun,F] = deformed_consensus(A)
%DEFORMED_CONSENSUS  n-by-n QEP from multi-agent systems theory.
%  [COEFFS,FUN,F] = nlevp('deformed_consensus',A) constructs an n-by-n 
%  quadratic matrix polynomial (the deformed Laplacian),
%  Q(lambda) = lambda^2*(D - eye(n)) - lambda*A + eye(n) arising in the study
%  of the stability properties of the deformed consensus protocol in control theory.
%  D and A are respectively the degree matrix and (0,1)-adjacency matrix of an
%  undirected graph with n vertices, describing the information flow among the
%  agents in the in the network. Since the graph is undirected, A is an n-by-n
%  real symmetric matrix.
%  The matrices are returned in a cell array: COEFFS = {eye(n), -A, D-eye(n)}.
%  FUN is a function handle to evaluate the monomials 1, lambda, lambda^2
%  and their derivatives.
%  F is the function handle eye(n) - lambda*A + lambda^2*(D-eye(n)).
%  This problem has the properties: pep, qep, real, symmetric, scalable,
%  parameter-dependent.

%  Reference:
%  F. Morbidi, The deformed consensus protocol, Automatica,
%  49(10): 3049-3055, 2013.
%
%  This function has been contributed by F. Morbidi,
%  Universite de Picardie Jules Verne, July 17, 2015.


if nargin < 1
   load deformed_consensus;
end
if ~isreal(A), error('Input matrix must be real.');end
if ~isequal(A,A'), error('Input matrix must be symmetric.');end

[m,n] = size(A);
if m ~= n
   disp('Error: the input matrix A is not square')
end

% Computation of the degree matrix D from the adjacency matrix A
D = A*ones(n,1);

% Coefficient matrices of the QEP
M = diag(D) - eye(n);
C = -A;
K = eye(n);

coeffs = {K, C, M};

fun = @(lam) nlevp_monomials(lam,2);
F = nlevp_handleQEP(coeffs);

end
