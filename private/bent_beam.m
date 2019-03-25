function [coeffs,fun,F] = bent_beam
%BENT_BEAM  6-by-6 NEP from a bent beam model.
%  [COEFFS, FUN, F] = nlevp('bent_beam') constructs a 6-by-6
%  nonlinear eigenvalue problem A(lam) = 
%  A11*cosh(alpha(lam)*l) + A12*alpha(lam)*cosh(alpha(lam)*l)
%   + A13*alpha(lam)^2*cosh(alpha(lam)*l)
%   + A21*cos(alpha(lam)*l) + A22*alpha(lam)*cos(alpha(lam)*l)
%   + A23*alpha(lam)^2*cos(alpha(lam)*l)
%   + A31*sinh(alpha(lam)*l) + A32*alpha(lam)*sinh(alpha(lam)*l)
%   + A33*alpha(lam)^2*sinh(alpha(lam)*l)
%   + A41*sin(alpha(lam)*l) + A42*alpha(lam)*sin(alpha(lam)*l)
%   + A43*alpha(lam)^2*sin(alpha(lam)*l)
%   + A51*cos(beta(lam)*l) + A52*beta(lam)*cos(beta(lam)*l)
%   + A61*sin(beta(lam)*l) + A62*beta(lam)*sin(beta(lam)*l),
%  where alpha(lam) = k1*lam^(1/2) and beta(lam) = k2*lam and
%  l, k1, k2 are specific physical constants.
%  The matrices are returned in a cell array: COEFFS =
%  {A11, A12, A13, A21, A22, A23, A31, A32, A33, A41, A42, A43,
%  A51, A52, A61, A62} 
%  FUN is a array function that evaluates the 16 functions related to
%  the matrices in COEFFS.
%  F returns the function handle that evaluates the NEP.

%  Reference:
%  N. K. Jain, K. Singhal, and K. Huseyin. On roots of functional lambda 
%  matrices. Comput. Methods Appl. Mech. Engrg.,40(3):277–292, 1983.

l = 23.5;
EI = 38.92e3;
m = 1.833e-4;
tau = 4.93e4;

% Building the coefficient matrices
coeffs = cell(1,16);

A11 = zeros(6,6);
A11(1,1) = 1; A11(1,3) = -1;
A11(3,2) = 1; A11(3,4) = 1;
coeffs{1} = A11;

A12 = zeros(6,6);
A12(2,4) = 1; A12(6,2) = 1;
coeffs{2} = A12;

A13 = zeros(6,6);
A13(4,3) = EI; A13(5,1) = EI;
coeffs{3} = A13;

A21 = A11;
A21(1,1) = -1; A21(3,4) = -1;
coeffs{4} = A21;

A22 = A12;
A22(6,2) = -1;
coeffs{5} = A22;

A23 = A13;
A23(4,3) = -EI;
coeffs{6} = A23;

A31 = zeros(6,6);
A31(1,2) = 1; A31(1,4) = -1;
A31(3,1) = 1; A31(3,3) = 1;
coeffs{7} = A31;

A32 = zeros(6,6);
A32(2,3) = 1; A32(6,1) = 1;
coeffs{8} = A32;

A33 = zeros(6,6);
A33(4,4) = EI; A33(5,2) = EI;
coeffs{9} = A33;

A41 = A31;
A41(1,2) = -1; A41(3,1) = -1;
coeffs{10} = A41;

A42 = A32;
A42(2,3) = -1;
coeffs{11} = A42;

A43 = A33;
A43(4,4) = -EI;
coeffs{12} = A43;

A51 = zeros(6,6);
A51(6,6) = -1;
coeffs{13} = A51;

A52 = zeros(6,6);
A52(4,5) = -tau;
coeffs{14} = A52;

A61 = zeros(6,6);
A61(2,5) = 1;
coeffs{15} = A61;

A62 = zeros(6,6);
A62(5,6) = tau;
coeffs{16} = A62;
% End of building coefficient matrices

fun = @bent_beam_fun;
%Building the function_handle F
m1EI4 = (m/EI)^(0.25);
lm1EI4 = l*m1EI4;
m1t2 = (m/tau)^(0.5);
lm1t2 = l*m1t2;
mEI2 = (m*EI)^0.5;
mt2 = (m*tau)^0.5;

F1 = @(z) [cosh(lm1EI4*sqrt(z))-cos(lm1EI4*sqrt(z)) ...
     sinh(lm1EI4*sqrt(z))-sin(lm1EI4*sqrt(z)) -(cosh(lm1EI4*sqrt(z))+cos(lm1EI4*sqrt(z)))...
	 -(sinh(lm1EI4*sqrt(z))+sin(lm1EI4*sqrt(z))) 0 0];

F2 = @(z)   [0 0 m1EI4*sqrt(z)*(sinh(lm1EI4*sqrt(z))-sin(lm1EI4*sqrt(z)))...
     m1EI4*sqrt(z)*(cosh(lm1EI4*sqrt(z))+cos(lm1EI4*sqrt(z))) sin(lm1t2*z) 0];
    
F3 = @(z) [sinh(lm1EI4*sqrt(z))-sin(lm1EI4*sqrt(z)) ...
     cosh(lm1EI4*sqrt(z))+cos(lm1EI4*sqrt(z))...
     sinh(lm1EI4*sqrt(z))+sin(lm1EI4*sqrt(z)) cosh(lm1EI4*sqrt(z))-cos(lm1EI4*sqrt(z)) 0 0];

F4 = @(z) [0 0 mEI2*z*(cosh(lm1EI4*sqrt(z))-cos(lm1EI4*sqrt(z)))...
     mEI2*z*(sinh(lm1EI4*sqrt(z))-sin(lm1EI4*sqrt(z))) -mt2*z*cos(lm1t2*z) 0];

F5 = @(z) [ mEI2*z*(cosh(lm1EI4*sqrt(z))+cos(lm1EI4*sqrt(z)))...
	 mEI2*z*(sinh(lm1EI4*sqrt(z))+sin(lm1EI4*sqrt(z))) 0 0 0 mt2*z*sin(lm1t2*z)];

F6 = @(z) [m1EI4*sqrt(z)*(sinh(lm1EI4*sqrt(z))+sin(lm1EI4*sqrt(z)))...
     m1EI4*sqrt(z)*(cosh(lm1EI4*sqrt(z))-cos(lm1EI4*sqrt(z))) 0 0 0 -cos(lm1t2*z)];

F = @(z) [F1(z); F2(z); F3(z); F4(z); F5(z); F6(z)];

end

function F = bent_beam_fun(z)
    
l= 23.5;
EI = 38.92e3;
m = 1.833e-4;
tau = 4.93e4;
      
% Building the functions
fun11 =  cosh(l*(m/EI)^(0.25)*sqrt(z));
fun12 =  (m/EI)^(0.25)*sqrt(z)*cosh(l*(m/EI)^(0.25)* ...
                                        sqrt(z));
fun13 =  (m/EI)^(0.5)*z*cosh(l*(m/EI)^(0.25)*sqrt(z));
fun21 =  cos(l*(m/EI)^(0.25)*sqrt(z));
fun22 =  (m/EI)^(0.25)*sqrt(z)*cos(l*(m/EI)^(0.25)* ...
                                       sqrt(z));
fun23 =  (m/EI)^(0.5)*z*cos(l*(m/EI)^(0.25)*sqrt(z));
fun31 =  sinh(l*(m/EI)^(0.25)*sqrt(z));
fun32 =  (m/EI)^(0.25)*sqrt(z)*sinh(l*(m/EI)^(0.25)* ...
                                        sqrt(z));
fun33 =  (m/EI)^(0.5)*z*sinh(l*(m/EI)^(0.25)*sqrt(z));
fun41 =  sin(l*(m/EI)^(0.25)*sqrt(z));
fun42 =  (m/EI)^(0.25)*sqrt(z)*sin(l*(m/EI)^(0.25)* ...
                                       sqrt(z));
fun43 =  (m/EI)^(0.5)*z*sin(l*(m/EI)^(0.25)*sqrt(z));
fun51 =  cos(l*(m/tau)^(0.5)*z);
fun52 =  (m/tau)^(0.5)*z*cos(l*(m/tau)^(0.5)*z);
fun61 =  sin(l*(m/tau)^(0.5)*z);
fun62 =  (m/tau)^(0.5)*z*sin(l*(m/tau)^(0.5)*z);
    
F =  [fun11, fun12, fun13, fun21, fun22, fun23, fun31, fun32, fun33,  ...
      fun41, fun42, fun43, fun51, fun52, fun61, fun62];
% End of building the functions
end
