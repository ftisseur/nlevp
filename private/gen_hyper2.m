function [coeffs,fun,F,sol] = gen_hyper2(n,E,V)
%GEN_HYPER2    Hyperbolic QEP constructed from prescribed eigenpairs.
%  [COEFFS,FUN,F,SOL] = nlevp('gen_hyper2',N,S) generates an N-by-N
%  random hyperbolic quadratic matrix polynomial
%  Q(lambda) = lambda^2*A + lambda*B + C
%  with eigenvalues and eigenvectors returned in SOL.EVAL and SOL.EVEC.
%  [COEFFS,FUN,F,SOL] = nlevp('gen_hyper2',E,S) generates Q(lambda)
%  with eigenvalues E, where E must be a vector of length 2*N of real
%  numbers ordered so that min(E(1:N)) > max(E(N+1:2*N)).
%  [COEFFS,FUN,F,SOL] = nlevp('gen_hyper2',E,V,S) is similar and produces
%  Q(lambda) with eigenvectors V, where V is an N-by-2*N matrix
%  satisfying V*[EYE(N) 0; 0 -EYE(N)]*V' = 0.
%  The optional parameter S (default: S = 0) is a nonnegative integer used to 
%  seed the random number generator. Set S = 'noseed' to leave the random 
%  number generator unseeded.
%  Defaults: N = 15, E = SORT(N*RANDN(2*N,1),'descend'), V = [V1 V1*U],
%  where V1 = RANDN(N) and U = ORTH(RANDN(N)).
%  The matrices are returned in a cell array: COEFFS = {C, B, A}.
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2
%  and their derivatives.
%  F is the function handle C + lambda*B + lambda^2*A.
%  This problem has the properties pep, qep, real, symmetric, hyperbolic,
%  parameter-dependent, scalable, solution, random.

%  References:
%  M. Al-Ammari and F. Tisseur, Hermitian matrix polynomials with real
%    eigenvalues of definite type. Part I: Classification. To appear in
%    Linear Algebra and Appl., 2011. doi:10.1016/j.laa.2010.08.035
%  C.-H. Guo, N. J. Higham, and F. Tisseur, Detecting and solving
%    hyperbolic quadratic eigenvalue problems. SIAM J. Matrix Anal.
%    Appl., 30(4):1593-1613, 2009.

if nargin < 1, n = 15; E = []; V = []; state = []; end
if nargin == 1
   if length(n) > 1
      E = n; n = length(E)/2;
   else
      E = [];
   end
   V = []; state = [];
end
if nargin == 2
  if length(n) > 1   % E, S or E, V case
    if length(E) > 1 % E, V case
       V = E; E = n; n = length(E)/2; state = 'noseed';
    else             % E, S case
       state = E; E = n; n = length(E)/2; V = [];
    end
  else               % N, S case
    state = E; E = []; V = [];
  end
end
if nargin == 3, state = V; V = E; E = n; n = length(E)/2; end
if n ~= round(n), error('E must have an even number of elements.'), end

if isempty(state)
  warning('NLEVP:random','Random numbers are now seeded by this function.')
  state = 0;
end

if ~strcmpi(state,'noseed')
   rng(state);
end

if ~isempty(V)
   if ~all(size(V) == [n 2*n])
      error('Third input argument has incorrect size.')
   end
   V1 = V(:,1:n);
   V2 = V(:,n+1:2*n);
   if norm(V1*V1'-V2*V2',1) > n*eps*(norm(V1,1)^2 + norm(V2,1)^2)
      error('Incorrect set of eigenvectors.')
   end
else
   V1 = orth(randn(n)); % Arbitrary, but orth. leads to better conditioned A.
   if nlevp_isoctave
      % Octave has no gallery function.
      [Q,R] = qr(randn(size(V1,2)));
      V2 = Q*diag(sign(diag(R)))*V1';
   else
      V2 = gallery('qmult',V1',1)';
   end
end

if ~isempty(E)
   E1 = E(1:n);
   E2 = E(n+1:2*n);
   if min(E1) < max(E2) || ~isreal(E)
      error('Incorrect set of eigenvalues.')
   end
else
   gap = 0;
   while gap < n*eps
     E = sort(n*randn(2*n,1));
     E2 = E(1:n); E1 = E(n+1:2*n);
     gap = abs(E1(n)-E2(1));
   end
end

Ainv = V1*diag(E1)*V1'-V2*diag(E2)*V2'; % Ainv is positive definite.
Ainv = (Ainv+Ainv')/2; % To ensure symmetry.
A = inv(Ainv);
B = -A*(V1*diag(E1.^2)*V1'-V2*diag(E2.^2)*V2')*A;
B = .5*(B+B');
C = -A*(V1*diag(E1.^3)*V1'-V2*diag(E2.^3)*V2')*A + B*Ainv*B;
C = .5*(C+C');

coeffs = {C, B, A};
fun = @(lam) nlevp_monomials(lam,2);
F = nlevp_handleQEP(coeffs);
sol.evec = [V1 V2];
sol.eval = [E1; E2];


end
