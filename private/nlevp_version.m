function varargout = nlevp_version
%NLEVP_VERSION   Version helper function for NLEVP.
%  See HELP NLEVP for syntax.

load nlevp_version.mat % Load structure v.

if nargout == 0
    fprintf(['This is the NLEVP collection of nonlinear eigenvalue ',...
             'problems version %s.\n',...
             'It was released %s and contains %d problems.\n'],...
             v.number, v.date, v.problemcount);
else
    varargout{1} = v;
end

end
