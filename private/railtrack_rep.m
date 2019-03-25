function [coeffs,fun,F,xcoeffs] = railtrack_rep
%RAILTRACK_REP  REP from study of vibration of rail tracks.
%  [COEFFS,FUN,F,XCOEFFS] = nlevp('railtrackREP') generates the coefficient
%  matrices of the railtrack problem, a quadratic eigenvalue problem
%      (lambda C +  B + A/lambda)*x = 0,
%  which models vibration of rail tracks under the excitation of high speed
%  trains.  This problem is T-palindromic, i.e., C = A.' and B = B.',
%  and has dimension 1005.
%  The matrices are sparse and returned in a cell array: COEFFS = {A, B, C}.
%  Notice that this is the rational formulation of nlevp('railtrack').
%  F is the function handle lambda C + B + \lambda^{-1}A.
%  FUN is a function handle to evaluate the monomials \lambda^{-1}, 1, lambda,
%  and their derivatives.
%  XCOEFFS is the cell {L1 1 L3; R1 B R3} that exploits the low rank of A =
%  L1*R1' and of C = L3*R3'.
%  This problem has the properties rep, sparse, low-rank.

%  References:
%  A. Hilliges. Numerische L\"osung von quadratischen Eigenwertproblemen mit 
%  Anwendung in der Schienendynamik. Master's thesis, TU Berlin, 2004.
%
%  Andreas Hilliges, Christian Mehl, and Volker Mehrmann. On the solution of 
%  palindromic eigenvalue problems. In P. Neittaanmaki, T. Rossi, S. Korotov,
%  E. Onate, J. Periaux, and D. Knorzer, editors, Proceedings of the European 
%  Congress on Computational Methods in Applied Sciences and Engineering 
%  (ECCOMAS 2004), Jyvaskyla, Finland, 2004.
%  http://www.mit.jyu.fi/eccomas2004/proceedings/proceed.html
%
%  D. Steven Mackey, Niloufer Mackey, Christian Mehl, and Volker Mehrmann.
%  Structured polynomial eigenvalue problems: Good vibrations from
%  good linearizations. SIAM J. Matrix Anal. Appl., 28(4):1029-1051, 2006.

load('railtrack.mat')
coeffs = {sA.', sB, sA};
fun = @(lam) railtrackREP_fun(lam);
F = @(lam) coeffs{1}/lam + coeffs{2}  + lam*coeffs{3};
xcoeffs = {R3, 1, L3; L3', sB, R3'};

end


function varargout = railtrackREP_fun(lam)

lam = lam(:);
n = length(lam);
varargout{1} = [lam.^(-1), ones(n,1), lam];
if nargout >= 2
    f = -(lam.^(-2));
    varargout{2} = [f, zeros(n,1),ones(n,1)];
    for i=2:nargout-1
        f = -i*f./lam;
        varargout{i+1} = [f, zeros(n,2)];
    end
end
end
