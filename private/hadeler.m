function [coeffs,fun,F] = hadeler(n,alpha)
%HADELER  NEP due to Hadeler.
%  [COEFFS,FUN,F] = nlevp('hadeler',N,ALPHA) returns matrices B, A_2,
%  and A_0 for the nonlinear eigenvalue problem of dimension N defined by
%  T(lam)*x = [(exp(lam)-1)*B + lam^2*A_2 - A_0]*x = 0.
%  A_0 = ALPHA*EYE(N), where the defaults are N=8, ALPHA = 100.
%  The matrices are returned in a cell array: COEFFS = {A_0,A_2,B}.
%  FUN is a function handle to evaluate the functions -1,lam^2,
%  and exp(lam)-1 and their derivatives.
%  F is the function handle T(lam). 
%  This problem has the properties nep, real, symmetric, scalable.

%  References:
%  K. P. Hadeler. Mehrparametrige und nichtlineare Eigenwertaufgaben.
%  Arch. Rational Mech. Anal., 27(4):306-328, 1967.

if nargin < 1 || isempty(n), n = 8; end
if nargin < 2 || isempty(alpha), alpha = 100; end

i = 1:n;
I = i(ones(n,1),:);

A_0 = alpha*eye(n);
A_2 = n*eye(n) + 1./(I+I');
B   = (n+1 - max(I,I')) .* (i'*i);

coeffs = {A_0,A_2,B};
fun = @hadeler_fun;
F = @(lam) (exp(lam)-1)*B + lam^2*A_2 - A_0;

end

function varargout = hadeler_fun(lam)
lam = lam(:);
n = length(lam);
explam = exp(lam);

varargout{1} = [-ones(n,1),lam.^2,explam-1];
if nargout >= 2
    varargout{2} = [zeros(n,1),2*lam,explam];
end
if nargout >= 3
    varargout{3} = [zeros(n,1),2*ones(n,1),explam];
end
for i=4:nargout
    varargout{i} = [zeros(n,2),explam];
end

end
