function [coeffs, fun, F] = neuron_dde(param)
%NEURON_DDE    2-by-2 NEP from a neural-network DDE.
% [coeffs, fun, F] = nlevp('neuron_dde', param) generates the 2-by-2
% nonlinear matrix function lam*EYE(2) + B0 - exp(-lam*taus)*A0 -
% exp(-lam*tau1)*A1 - exp(-lam*tau2)*A2, which is the characteris equation
% of a delay-differential equation from a neural network with two neurons.
% COEFFS is the cell of the matrix coefficients {EYE(2), B0, A0, A1, A2},
% where B0 is k*eye(2), A0 is beta*eye(2), A1 is [0 a12; 0 0] and A2 is
% [0 0; a21 0]. The optional parameter 'param' should contain the values
% [k, beta, a12, a21, taus, tau1, tau2]. The default behaviour is
% [0.5,-1,1,2.5,0.01,1,1]. 
% FUN is a function handle to evaluate the functions [lam, 1,
% -exp(-lam*taus), -exp(-lam*tau1), -exp(-lam*tau2)] and their derivatives.
% F is a function handle to evaluate the matrix-valued function lam*EYE(2)
% + B0 - exp(-lam*taus)*A0 - exp(-lam*tau1)*A1 - exp(-lam*tau2)*A2.
% The problem has the properties nep, real, parameter-dependent

% Reference: Section 2 of L. P. Shayer, and S. A. Campbell,
% "Stability, bifurcation and multistability in a system of two coupled
% neurons with multiple time delays", SIAM J. APPL. MATH.  Vol. 61, No. 2,
% pp. 673–700, 2000. DOI: https://doi.org/10.1137/S0036139998344015

if nargin < 1 || isempty(param)
param = [0.5,-1,1,2.5,0.01,1,1];
end
% setting all the parameters
k = param(1);
beta = param(2);
a12 = param(3);
a21 = param(4);
taus = param(5);
tau1 = param(6);
tau2 = param(7);

I = eye(2);
coeffs = cell(1,5);
coeffs{1} = I;
coeffs{2} = k*I;
coeffs{3} = beta*I;
coeffs{4} = [0 1 ; 0 0];
coeffs{5} = [0 0 ; 1 0];

F = @(lam) lam*coeffs{1} + coeffs{2} -exp(-lam*taus)*coeffs{3} ...
    - exp(-lam*tau1)*coeffs{4} - exp(-lam*tau2)*coeffs{5};
fun = @(lam) neuron_dde_fun(lam, taus, tau1, tau2);
end


function varargout = neuron_dde_fun(lam, taus, tau1, tau2)

tau = [taus, tau1, tau2];
lam = lam(:);
l = length(lam);

varargout{1} = [lam ones(l, 1) -exp(-lam*tau)];
if nargout > 1
   deriv = exp(-lam*tau)*diag(tau);
   varargout{2} = [ones(l,1)  zeros(l,1) deriv];
   for j = 3:nargout
      deriv = -deriv*diag(tau);
       varargout{j} = [zeros(l,2) deriv];
   end
end
end