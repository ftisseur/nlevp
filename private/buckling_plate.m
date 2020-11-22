function [coeffs,fun,F] = buckling_plate
%BUCKLING_PLATE  3-by-3 NEP from a buckling plate model.
%  [COEFFS,FUN,F] = nlevp('buckling_plate') generates a 3-by-3 nonlinear matrix 
%  function  A_0 + f*A_1 + g*A_2, where f and g are ratios of trigonometric 
%  functions that model the buckling of a plate.
%  The problem is highly irregular: it is not defined in 0 and in many other points.
%  The matrices are returned in a cell array: COEFFS = {A_0, A_1, A_2}.
%  FUN is a function handle to evaluate the functions 1, f and g.
%  F is the function handle that evaluates the NEP.
%  This problem has the properties nep, real, symmetric.

%  This problem has been sent us by Melina Freitag, University of Bath.

A_0 = [10 0 2; 0 4 2; 2 2 8];
A_1 = [1  0 0; 0 1 0; 0 0 0];
A_2 = [0  1 0; 1 0 0; 0 0 0];

coeffs = {A_0, A_1, A_2};

fun = @(lam) buckling_plate_fun(lam);

%Building F
f = @(lam) lam.*(1 - 2*lam.*cot(2*lam))./(tan(lam) - lam);
g = @(lam) (lam.*(2*lam - sin(2*lam)))./(sin(2*lam).*(tan(lam) - lam));
F = @(lam) [f(lam)+10 g(lam) 2; g(lam) f(lam)+4 2; 2 2 8];

end


function F = buckling_plate_fun(lam)

lam = lam(:);
n = length(lam);
f = lam.*(1 - 2*lam.*cot(2*lam))./(tan(lam) - lam);
g =(lam.*(2*lam - sin(2*lam)))./(sin(2*lam).*(tan(lam) - lam));
F = [ones(n,1) f g];

end
