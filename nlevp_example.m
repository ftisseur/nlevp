function nlevp_example(fname)
%NLEVP_EXAMPLE  Run POLYEIG on PEP problems from NLEVP.
% NLEVP_EXAMPLE solves all the not-too-large PEP problems in NLEVP
% by POLYEIG, sending output to the screen.
% NLEVP_EXAMPLE(fname) directs partial output to the file named fname
% (intended for generating output for NLEVP paper).

if nargin == 0
   fid = 1;
else
   fid = fopen(fname,'w');
end
s_rand = warning('off', 'NLEVP:random');  % For gen_hyper2.

nmax = 500;
probs = nlevp('query','pep');
nprobs = length(probs);
nprobs_total = length(nlevp('query','problems'));
fprintf(fid,'NLEVP contains %2.0f problems in total,\n', nprobs_total);
fprintf(fid,'of which %2.0f are polynomial eigenvalue problems (PEPs).\n', nprobs);
fprintf(fid,'Run POLYEIG on the PEP problems of dimension at most %2.0f:\n\n',nmax);

fprintf(fid,'             Problem   Dim  Max and min magnitude of eigenvalues\n');
fprintf(fid,'             -------   ---  ------------------------------------\n');
m2 = 7;
m1 = ceil(nprobs/m2);
j = 1;
for i=1:nprobs
   if fid ~= 1 && i == 9
      fprintf(fid,'                 ...\n');
      fid_save = fid;
      fid = 1;  % Omit output from this point on when writing to file.
   end
   coeffs = nlevp(probs{i});
   [n, nc] = size(coeffs{1});
   if n >= nmax
      fprintf(fid,'%20s   %3.0f is a PEP but is too large for this test.\n', ...
              probs{i}, n);
   elseif n ~= nc
      fprintf(fid,'%20s   %3.0f is a PEP but is nonsquare.\n', probs{i}, n);
   else

      % POLYEIG will convert sparse input matrices to full.
      e = polyeig(coeffs{:});
      fprintf(fid,'%20s   %3.0f  %9.2e, %9.2e\n', ...
               probs{i}, n, max(abs(e)), min(abs(e)));
      subplot(m1,m2,j)
      plot(real(e), imag(e),'.')
      title(probs{i},'Interpreter','none')
      % Tweaks.
      if strcmp(probs{i},'sign1'), ylim([-1 1]*1.5), end
      if strcmp(probs{i},'damped_beam')
         title(['         ' probs{i}],'Interpreter','none')
      end
      if strcmp(probs{i},'relative_pose_6pt')
         title(['          ' probs{i}],'Interpreter','none')
      end
      if strcmp(probs{i},'speaker_box') || strcmp(probs{i},'intersection')
         title(['     ' probs{i}],'Interpreter','none')
      end
      j = j+1;
   end
end

if nargin > 0, fclose(fid_save); end
warning(s_rand)

end
