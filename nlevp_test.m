function nlevp_test
%NLEVP_TEST    Test properties of NLEVP collection.
%  For use by developers, but can also be used as a test for correctness
%  of installation.

% This needs to be a function so that NLEVP_ISOCTAVE can be called.

disp('Testing the NLEVP collection')
if nlevp_isoctave
   disp('Running Octave.')
else   
   disp('Running MATLAB.')
end
numerror = 0;
% Warnings are not currently needed, but might be in future if we need
% to test a property for which we have only necessary conditions.
% numwarning = 0;
try
    peps = nlevp('query','pep');
    probs = nlevp('query','problems');
catch % exception
      % "exception" commented out for Octave compatibility.
    disp('NLEVP does not work!!')
    return
end

s_rand = warning('off', 'NLEVP:random');  % For gen_hyper2.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Testing generation of all problems')
probsNew = probs;
for i=1:length(probs)
    try
        [coeffs,fun,~] = nlevp(probs{i});
        % Try to evaluate T(2).
        f = fun(2);
        T = coeffs{1}*f(1);
        for k=2:length(coeffs)
            T = T + coeffs{k}*f(k);
        end
    catch % exception
        warning('NLEVP:warning_no_generation',...
                'Error generating problem %s.', probs{i});
        probsNew = setdiff(probsNew,probs{i});
        peps = setdiff(peps,probs{i});
        numerror = numerror + 1;
    end
end
probs = probsNew;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Testing fun and F returning the same problem')
for i=1:length(probs)
  [coeffs,fun,F] = nlevp(probs{i});
  f = fun(2.3);
  T = coeffs{1}*f(1);
  for k=2:length(coeffs)
    T = T + coeffs{k}*f(k);
  end
  thresh = 1e-11;
  if norm(F(2.3)-T, 1)/norm(F(2.3), 1) > thresh
    warning('NLEVP:warning_diff_problems',...
            'Problem %s returns different nlevp in fun and F.', probs{i});
    numerror = numerror + 1;
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Testing T-palindromicity')
probs = nlevp('query','T-palindromic');
probs = intersect(probs,peps);
for i=1:length(probs)
    coeffs = nlevp(probs{i});
    k = length(coeffs);
    for j=1:ceil(k/2)
        if ~isequal(coeffs{j},coeffs{k-j+1}.')
            warning('NLEVP:warning_not_palindromic',...
                'Problem %s is claimed to be T-palindromic, but is not.',...
                probs{i})
            numerror = numerror + 1;
            break;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Testing *-palindromicity')
probs = nlevp('query','*-palindromic');
probs = intersect(probs,peps);
for i=1:length(probs)
    coeffs = nlevp(probs{i});
    k = length(coeffs);
    for j=1:ceil(k/2)
        if ~isequal(coeffs{j},coeffs{k-j+1}')
            warning('NLEVP:warning_not_palindromic',...
                'Problem %s is claimed to be *-palindromic, but is not',...
                probs{i})
            numerror = numerror + 1;
            break;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Testing T-anti-palindromicity')
probs = nlevp('query','T-anti-palindromic');
probs = intersect(probs,peps);
for i=1:length(probs)
    coeffs = nlevp(probs{i});
    k = length(coeffs);
    for j=1:ceil(k/2)
        if ~isequal(coeffs{j},-coeffs{k-j+1}.')
            warning('NLEVP:warning_not_palindromic',...
                ['Problem %s is claimed to be T-anti-palindromic, ',...
                'but is not.'], probs{i})
            numerror = numerror + 1;
            break;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Testing *-anti-palindromicity')
probs = nlevp('query','*-anti-palindromic');
probs = intersect(probs,peps);
for i=1:length(probs)
    coeffs = nlevp(probs{i});
    k = length(coeffs);
    for j=1:ceil(k/2)
        if ~isequal(coeffs{j},-coeffs{k-j+1}')
            warning('NLEVP:warning_not_palindromic',...
                ['Problem %s is claimed to be *-anti-palindromic, ',...
                'but is not.'], probs{i})
            numerror = numerror + 1;
            break;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Testing T-eveness')
probs = nlevp('query','T-even');
probs = intersect(probs,peps);
for i=1:length(probs)
    coeffs = nlevp(probs{i});
    k = length(coeffs);
    for j=1:k
        if ~isequal(coeffs{j},(-1)^(j-1)*coeffs{j}.')
            warning('NLEVP:warning_not_even',...
                'Problem %s is claimed to be T-even, but is not.',...
                probs{i})
            numerror = numerror + 1;
            break;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Testing *-eveness')
probs = nlevp('query','*-even');
probs = intersect(probs,peps);
for i=1:length(probs)
    coeffs = nlevp(probs{i});
    k = length(coeffs);
    for j=1:k
        if ~isequal(coeffs{j},(-1)^(j-1)*coeffs{j}')
            warning('NLEVP:warning_not_even',...
                'Problem %s is claimed to be *-even, but is not.',...
                probs{i})
            numerror = numerror + 1;
            break;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Testing T-oddness')
probs = nlevp('query','T-odd');
probs = intersect(probs,peps);
for i=1:length(probs)
    coeffs = nlevp(probs{i});
    k = length(coeffs);
    for j=1:k
        if ~isequal(coeffs{j},(-1)^(j)*coeffs{j}.')
            warning('NLEVP:warning_not_odd',...
                'Problem %s is claimed to be T-odd, but is not.',...
                probs{i})
            numerror = numerror + 1;
            break;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Testing *-oddness')
probs = nlevp('query','*-odd');
probs = intersect(probs,peps);
for i=1:length(probs)
    coeffs = nlevp(probs{i});
    k = length(coeffs);
    for j=1:k
        if ~isequal(coeffs{j},(-1)^(j)*coeffs{j}')
            warning('NLEVP:warning_not_odd',...
                'Problem %s is claimed to be *-odd, but is not.',...
                probs{i})
            numerror = numerror + 1;
            break;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Testing sparsity')
probs = nlevp('query','sparse');
for i=1:length(probs)
    coeffs = nlevp(probs{i});
    for j = 1:length(coeffs)
        if ~issparse(coeffs{j})
            warning('NLEVP:warning_not_sparse',...
                    'Problem %s is claimed to be sparse ', ...
                    'but coefficient matrix %g is not.',...
                    probs{i}, j)
            numerror = numerror + 1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Testing random problems')
% Test that problem is fixed for default parameters.
probs = nlevp('query','random');
for i=1:length(probs)
    coeffs1 = nlevp(probs{i});
    coeffs2 = nlevp(probs{i});
    if ~isequal(coeffs1,coeffs2)
            warning('NLEVP:warning_not_fixed',...
                    ['Problem %s returns different matrices on different', ...
                    ' calls.'], probs{i})
            numerror = numerror + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Testing symmetry')
probs = nlevp('query','symmetric');
for i=1:length(probs)
    T = nlevp('eval',probs{i},rand(1));
    if ~isequal(T,T.')
        warning('NLEVP:warning_not_symmetric',...
            'Problem %s is claimed to be symmetric, but is not.',...
            probs{i})
        numerror = numerror + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Testing hyperbolicity')
probs = nlevp('query','hyperbolic');
for i=1:length(probs)
    coeffs = nlevp(probs{i});

    % This code not now needed because it is only a necessary test compared
    % with the nec. and suff. test below.  But it is valid for any degree.
    x = randn(length(coeffs{1}),1);
    if x'*coeffs{2}*x <= (x'*coeffs{1}*x) * (x'*coeffs{3}*x) ...
       || min(eig(coeffs{3})) < 0
        warning('NLEVP:warning_not_hyperbolic',...
            'Problem %s is claimed to be hyperbolic, but is not.',...
            probs{i})
        numerror = numerror + 1;

    % This code assumes a quadratic - can rewrite for general degree
    % as and when necessary. \cite[Thm.~3.4, P1]{alti10}.
    if length(coeffs) ~= 3
        warning('NLEVP:warning_hyperbolic_test_invalid',...
                'We neede to extend the hyperbolic test for higher degrees!')
    end
    [X,e] = polyeig(coeffs{:});
    [e,ind] = sort(e); X = X(:,ind);
    R = 2*coeffs{3}*X*diag(e)+coeffs{2}*X;
    nn = length(e);
    s = zeros(nn,1);
    for j = 1:nn, s(j) = sign(X(:,j)'*R(:,j)); end
    test2 = all(s(1:nn/2) == -ones(nn/1,1));
    test3 = all(s(nn/2+1:nn) == ones(nn/2,1));
    if any(imag(e(:))) || abs(e(nn/2)-e(nn/2+1))<eps || ~test2 || ~test3
        warning('NLEVP:warning_not_hyperbolic',...
            'Problem %s is claimed to be hyperbolic, but is not.',...
             probs{i})
        numerror = numerror + 1;
    end

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Testing Hermitian structure')
probs = nlevp('query','hermitian');
for i=1:length(probs)
    T = nlevp('eval',probs{i},rand(1));
    if ~isequal(T,T')
        warning('NLEVP:warning_not_hermitian',...
            'Problem %s is claimed to be hermitian, but is not.',...
            probs{i})
        numerror = numerror + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Testing realness')
probs = nlevp('query','real');
for i=1:length(probs)
    T = nlevp('eval',probs{i},rand(1));
    if ~all(all(isreal(T)))
        warning('NLEVP:warning_not_real',...
            'Problem %s is claimed to be real, but is not.', probs{i})
        numerror = numerror + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Testing nonsquareness')
probs = nlevp('query','nonsquare');
for i=1:length(probs)
    coeffs = nlevp(probs{i});
    k = length(coeffs);
    for j=1:k
        if diff(size(coeffs{j})) == 0
            warning('NLEVP:warning_not_nonsquare',...
                'Problem %s is claimed to be nonsquare, but is not.',...
                probs{i})
            numerror = numerror + 1;
            break;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Testing squareness')
probs = nlevp('query','nonsquare');
probs = setdiff(nlevp('query','problems'),probs);
for i=1:length(probs)
    coeffs = nlevp(probs{i});
    k = length(coeffs);
    for j=1:k
        if diff(size(coeffs{j})) ~= 0
            warning('NLEVP:warning_not_square',...
                'Problem %s is claimed to be square, but is not.',...
                probs{i})
            numerror = numerror + 1;
            break;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Testing proportionally damping')
% Test sufficient condition for this property.
probs = nlevp('query','proportionally-damped');
for i=1:length(probs)
    coeffs = nlevp(probs{i});
    k = length(coeffs);
    if k ~= 3
        warning('NLEVP:warning_not_proportionally_damped',...
            ['Problem %s is claimed to be proportionally damped, ',...
            'but is not a qep.'], probs{i})
        numerror = numerror + 1;
        break;
    end

%     Essentially contained in the test below.
%     A = [coeffs{1}(:)  coeffs{3}(:) coeffs{2}(:)];
%     R = qr(A); R = triu(R(1:3,1:3));
%     if abs(R(3,3)) > eps(10*norm(full(R(:,3))))
%         warning('NLEVP:warning_not_proportionally_damped',...
%             ['Problem %s may not be proportionally damped, ',...
%             'because damping matrix is not linear combination ',...
%             'of the mass and stiffness matrices.'], probs{i})
%         numwarning = numwarning + 1;
%     end

    % A sufficent condition.  \cite[Thm.~2]{laza09}.
    % Requires nonsingular coeffs{3}, which is always the case so far.
    n = length(coeffs{1});
    A = (coeffs{2}/coeffs{3})*coeffs{1}-(coeffs{1}/coeffs{3})*coeffs{2};
    if norm(A,1) > n*eps
        warning('NLEVP:warning_not_proportionally_damped',...
            ['Problem %s may not be proportionally damped as ',...
            'it does not satisfy the Caughey and O''Kelly condition.'],...
            probs{i})
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Testing given solutions')
probs = nlevp('query','solution');
% We need two lists of problems: the ones with the optional xcoeffs and the
% ones without.
probsLowRank = nlevp('query','low-rank');
probs1 = intersect(probs, probsLowRank);
probs2 = setdiff(probs, probsLowRank);

for i=1:length(probs1)
    [coeffs,fun,~,~,sol] = nlevp(probs1{i});
    k = length(coeffs);
    n = size(coeffs{1},1);
    m = length(sol.eval);
    R = zeros(n,m);
    W = zeros(n,m);
    lam = sol.eval; lam = lam(:);
    F = fun(lam);
    for j=1:k
        R = R + coeffs{j}*sol.evec*diag(F(:,j));
        W = W + abs(coeffs{j})*abs(sol.evec)*diag(abs(F(:,j)));
    end
    R = abs(R)./W; R(isnan(R)) = 0;
    if max(max(abs(R))) > 10^4*eps
        warning('NLEVP:warning_not_solution',...
                 ['Given solution of Problem %s ',...
                 'does not seem to be a solution.'], probs1{i})
        numerror = numerror + 1;
    end

end

for i=1:length(probs2)
    [coeffs,fun,~,sol] = nlevp(probs2{i});
    k = length(coeffs);
    n = size(coeffs{1},1);
    m = length(sol.eval);
    R = zeros(n,m);
    W = zeros(n,m);
    lam = sol.eval; lam = lam(:);
    F = fun(lam);
    for j=1:k
        R = R + coeffs{j}*sol.evec*diag(F(:,j));
        W = W + abs(coeffs{j})*abs(sol.evec)*diag(abs(F(:,j)));
    end
    R = abs(R)./W; R(isnan(R)) = 0;
    if max(max(abs(R))) > 10^4*eps
        warning('NLEVP:warning_not_solution',...
                 ['Given solution of Problem %s ',...
                 'does not seem to be a solution.'], probs2{i})
        numerror = numerror + 1;
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Testing low-rank')
probs = nlevp('query','low-rank');
for i = 1:length(probs)
    [coeffs, ~, ~, xcoeffs] = nlevp(probs{i});
    for j = 1:length(coeffs)
       thresh = norm(coeffs{j} - xcoeffs{1,j}*xcoeffs{2,j},inf)...
           /norm(coeffs{j},inf);
       if thresh > 1e-13
       warning(['NLEVP:warning_not_low-rank ',...
                   'Problem %s is claimed to be low-rank ', ...
                    'but coeffs %g is not equal to ', ...
                    'xcoeffs{1,%g}*xcoeffs{2,%g}'], probs{i}, j, j, j)
           numerror = numerror + 1;
       end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Testing tridiagonal')
probs = nlevp('query','tridiagonal');
for i=1:length(probs)
    coeffs = nlevp(probs{i});
    for j = 1:length(coeffs)
        if ~isbanded(coeffs{j},1,1)
            warning('NLEVP:warning_not_tridiagonal',...
                    'Problem %s is claimed to be tridiagonal ', ...
                    'but coefficient matrix %g is not.',...
                    probs{i}, j)
            numerror = numerror + 1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Testing banded')
probs = nlevp('query','banded');
lens = cellfun(@numel,probs);
pdlens = max(lens) - lens;
for i=1:length(probs)
    coeffs = nlevp(probs{i});
    N = size(coeffs{1},1);
    % we say that a matrix is banded
    % if the bandwidth is at most 10% of the size.
    max_bandwidth = N/5;
    bandwidth = 0;
    for j = 1:length(coeffs)
        k = 0;
        banded = 0;
        while ~banded && k <= max_bandwidth
            banded = isbanded(coeffs{j},k,k);
            k = k + 1;
        end
        bandwidth = max(k, bandwidth);
        if ~banded
            warning(['NLEVP:warning_not_banded: ',...
                'Problem %s is claimed to be banded ', ...
                'but coefficient matrix %g is not.'],...
                probs{i}, j)
            numerror = numerror + 1;
            continue
        elseif j == length(coeffs)
           fprintf('   Problem %s%s : bandwidth = %d%s size = %d\n',...
               probs{i}, repmat(' ',1,pdlens(i)), bandwidth,...
               repmat(' ',1,5-ceil(log10(bandwidth))), N) 
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('NLEVP collection tests completed.')
fprintf('***  Errors: %g\n', numerror)
% fprintf('***  Warnings: %g\n', numwarning)
warning(s_rand)

end
