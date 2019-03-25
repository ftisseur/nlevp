function [coeffs,fun,F] = bcc_traffic(n,s)
%BCC_TRAFFIC  QEP from stability analysis of chain of non-identical cars.
%  [COEFFS,FUN,F] = nlevp('bcc_traffic',KV,KD) generates an N-by-N monic
%  tridiagonal quadratic matrix polynomial Q(lambda) = lambda^2*I+ lambda*D + K,
%  where n >= 2, that arises from the stability analysis of a chain of N
%  non-identical vehicles under bilateral cruise control (BCC).
%  KV and KD are vectors of length N with positive entries corresponding to
%  the proportional gains and derivative gains of the cars.
%  [COEFFS,FUN,F] = nlevp('bcc_traffic',N,S) is similar and produces Q(lambda)
%  with randomly generated proportional gains and derivative gains.
%  kd(ind1) = .05+.1*randn(n1,1), kd(ind2) = .2+.2*rand(n2,1),
%  kv(ind2) = .05+.1*randn(n1,1), kv(ind1) = .2+.2*rand(n2,1), where
%  ind1 = ind(1:n1), ind2 = ind(n1+1:n), ind = randperm(n), n1+n2=N.
%  The default values are N = 32 (n1 = floor(N/2)).
%  The optional parameter S (default: S = 0) is a nonnegative integer used to 
%  seed the random number generator. Set S = 'noseed' to leave the random 
%  number generator unseeded.
%  The matrices are returned in a cell array: coeffs = {K, D, I}.
%  FUN is a function handle to evaluate the monomials 1,lambda,lambda^2
%  and their derivatives.
%  This problem has the properties pep, qep, real, parameter_dependent,
%  scalable, sparse, tridiagonal, banded.

%  References:
%  L. Wang, F. Tisseur, G. Strang, and B. K. P. Horn. Stability analysis of
%  a chain of non-identical vehicles under bilateral cruise control.
%  MIMS EPrint 2019.3, Manchester Institute for Mathematical Sciences,
%  The University of Manchester, UK, March 2019.

if nargin < 1 || isempty(n), n = 32; state = []; kv = []; kd = []; end
if nargin == 1
   if length(n) > 1
      error('n must be an integer')
   else
      state = []; kv = []; kd = [];
   end
end
if n <= 1, error('The integer n must be at least 2.'), end
if nargin == 2
   if length(n) > 1 % KV, KD case
      kv = n; kd = s; n = length(kv);
      if ~length(kd)==n, error('KV and KD must have same length'),end
    else             % N, S case
      state = s; kv = []; kd = [];
   end
end

if isempty(kv) || isempty(kd) %Construct kd and or kv

   if isempty(state)
     warning('NLEVP:random','Random numbers are now seeded by this function.')
     state = 0;
   end
   if ~strcmpi(state,'noseed')
     rng(state);
   end
   n1 = floor(n/2);
   %n2 = n-n1;
   index_car = randperm(n);
   index_1 = index_car(1: n1) ;
   index_2 = index_car(n1+1:n);
   if isempty(kv)
      kv = zeros(n,1);
      kv(index_1)= 0.3 + 0.2*rand(size(index_1));
      kv(index_2)= 0.1 + 0.1*(rand(size(index_2)) - 0.5);
   end;
   if isempty(kd)
      kd = zeros(n,1);
      kd(index_1)= 0.1 + 0.1*(rand(size(index_1)) - 0.5);
      kd(index_2)= 0.3 + 0.2*rand(size(index_2));
   end
end
if length(kd) ~= length(kv), error('kv and kd must have same length');end
if sum(kv > 0) ~= n, error('kv must be positive entries');end
if sum(kd > 0) ~= n, error('kd must be positive entries');end

S = diag(sparse(-ones(n-1,1)),-1)+ diag(sparse([2*ones(n-1,1);1]))+...
    diag(sparse(-ones(n-1,1)),1);

I = diag(sparse(ones(n,1)));
D = diag(sparse(kv))*S; %D = kv.*S;
K = diag(sparse(kd))*S; %K = kd.*S;
coeffs = {K, D, I};
fun = @(lam) nlevp_monomials(lam,2);
F =  nlevp_handleQEP(coeffs);

end
