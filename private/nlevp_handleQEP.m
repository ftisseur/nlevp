function F = nlevp_handleQEP(coeffs)
%NLEVP_HANDLE_QEP  builds a function_handle for QEPs
% F= nlevp_handle_QEP(K,D,M) builds the function_handle
% F = @(z) K + lam.*D + lam.^2.*M, where K, D, M are square
% matrices of the quadratic eigenvalue problem (QEP)
% K + lam*D + lam^2*M.
%
% This is a helper function. It is called by all the QEP problems. It
% should not be called directly by the user.

  F =  @(lam) coeffs{1} + lam.*(coeffs{2} + lam.*coeffs{3});

end
