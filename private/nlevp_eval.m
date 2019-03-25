function varargout = nlevp_eval(name,lam,varargin)
%NLEVP_EVAL     Evaluate nonlinear matrix function.
%  [T,TP,TPP,...] = nlevp('eval',NAME,LAMBDA,ARG1,ARG2,...)
%  evaluates the matrix function T and its derivatives TP, TPP,...
%  of problem NAME at the scalar LAMBDA.
%  ARG1, ARG2, ... are optional parameters depending on the problem.

if length(lam) > 1, error('Input argument LAMBDA must be a scalar.'), end

[coeffs,fun] = nlevp(name,varargin{:});
[F{1:nargout}] = fun(lam);

varargout = cell(1,nargout);
for k=1:nargout
    T = coeffs{1}*F{k}(1);
    for i=2:length(coeffs)
        T = T + coeffs{i}*F{k}(i);
    end
    varargout{k} = T;
end

end