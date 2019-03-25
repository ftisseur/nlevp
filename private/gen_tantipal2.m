function [coeffs,fun,F] = gen_tantipal2(n,T)
%GEN_TANTIPAL2 T-anti-palindromic QEP with eigenvalues on the unit circle.
%  [COEFFS,FUN,F] = nlevp('gen_tantipal2',N,S) generates an N-by-N
%  real T-anti-palindromic quadratic matrix polynomial
%  Q(lambda) = lambda^2*A + lambda*B - A', B = -B'
%  with eigenvalues randomly distributed on the unit circle.
%  [COEFFS,FUN] = nlevp('gen_tantipal2',T,S) generates Q(lambda)
%  with eigenvalues [cos(T)+i*sin(T);cos(T)-i*sin(T)], where
%  T is a vector of length N of distinct real numbers in (0,pi).
%  The optional parameter S (default: S = 0) is a nonnegative integer used 
%  to seed the random number generator (see HELP RAND).
%  Set S = 'noseed' to leave the random number generator unseeded.
%  Defaults: N = 16, T = a + (b-a).*rand(N,1), where a = 1e-4, b = pi - 1e-4.
%  The matrices are returned in a cell array: COEFFS = {-A', B, A}.
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2
%  and their derivatives.
%  F is the function handle -A' + lambda*B + lambda^2*A.
%  This problem has the properties pep, qep, real, T-anti-palindromic,
%  parameter-dependent, scalable, random.

%  References:
%  M. Al-Ammari, Analysis of structured polynomial eigenvalue problems.
%  PhD thesis, The University of Manchester, MIMS EPrint 2011.89,
%  Manchester Institute for Mathematical Sciences,
%  The University of Manchester, UK, 2011.

if nargin < 1, n = 16; T = []; state = []; end
if nargin == 1
   if length(n) > 1
      T = n; n = length(T);
   else
      T = [];
   end
   state = [];
end
if nargin == 2
  if length(n) > 1
        state = T; T = n; n = length(T);
  else
        state = T; T = [];
  end
end

if isempty(state)
  warning('NLEVP:random','Random numbers are now seeded by this function.')
  state = 0;
end

if ~strcmpi(state,'noseed')
   rng(state);
end

if isempty(T)
  tol = 1e-4; 
  a = tol; b = pi - tol;    % Exclude 0 and pi.
  T = a + (b-a).*rand(n,1); % Random angles.
end
if any(T <= 0) || any(T >= pi) 
   error('Elements of vector T must lie in open interval (0,1).')
end   
if any(T <= 1e-6) || any(T >= pi-1e-6) 
   warning('NLEVP:illcond',['Elements of T are within 1e-6 of 0 ' ...
            'or 1:ill conditioning may adversely effect the construction.'])
end   
e = cos(T) + 1i*sin(T); % Eigenvalues on upper part of unit circle.

% Map eigenvalues to the imaginary axis via Cayley transform and solve
% an inverse real T-even quadratic eigenvalue problem
% via the construction of a  a real Jordan triple (X,J,SX').

beta = imag((e-1)./(e+1));
b = zeros(2*n-1,1);
b(1:2:2*n-1) = beta;
J = diag(b,1)-diag(b,-1);   % Jordan matrix.
S = blkdiag(eye(n),-eye(n)); 
X1 = orth(rand(n));
if nlevp_isoctave
   % Octave has no gallery function.
   [Q,R] = qr(randn(n));
   X2 = Q*diag(sign(diag(R)))*X1';
else
   X2 = gallery('qmult',X1',1)';
end
X = [X1 X2];

% Now build the coefficient matrices of a real T-even quadratic.
JS = J*S;
A2 = inv(X*JS*X'); A2 = .5*(A2-A2');
JS = J*JS; temp = X*JS*X';
A1 = -A2*temp*A2; A1 = .5*(A1+A1');
A0 = -A2*(temp*A1+(X*(J*JS)*X')*A2) ; A0 = .5*(A0-A0');

% Apply inverse Cayley transform to get a real T-palindromic quadratic.
A = A2 + A1 + A0;
B = -2*A2 + 2*A0;

coeffs = {-A',.5*(B-B'),A};
fun = @(lam) nlevp_monomials(lam,2);
F = nlevp_handleQEP(coeffs);

end
