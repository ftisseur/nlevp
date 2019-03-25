function varargout = nlevp_monomials(lambda,k)
%NLEVP_MONOMIALS    Evaluate monomials and their derivatives.
%  [F,FP,FPP,...] = nlevp_monomials(LAM,K), for a column vector LAM,
%   generates
%    F = [1, LAM, LAM.^2, ..., LAM.^(K-1), LAM.^K],
%    FP = [0, 1, 2*LAM, 3*LAM.^2, ..., K*LAM.^(K-1)]'  (first derivative)
%  and higher derivatives FPP, ...
%  Thus for scalar LAM: F, FP, FPP, ... are row vectors of length K+1,
%  while for vectors LAM: F, FP, FPP, ... contain a row per element of LAM.
%
%  Define a polynomial P(x) = sum_{i=0:k} c(i)*x^i.
%  Then P(LAM) can be evaluated as P(LAM) = F*c.
%  The second and third derivatives can be evaluated as
%    P'(LAM) = FP*c, P''(LAM) = FPP*c,
%  and so on.

lambda = lambda(:);
n = length(lambda);

f = [ones(n,1),cumprod(repmat(lambda,1,k),2)];

varargout = cell(1,nargout);
varargout{1} = f;
for i=2:nargout
    f = [zeros(n,1),f(:,1:k).*repmat(1:k,n,1)];
    varargout{i} = f;
end
end