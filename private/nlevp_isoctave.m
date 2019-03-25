function k = nlevp_isoctave
%NLEVP_ISOCTAVE    Test whether Octave or MATLAB is running.
%  NLEVP_ISOCTAVE is TRUE if Octave is running and FALSE if MATLAB
%  is running.

k = false;
v = ver;
for j=1:length(v)
   if strfind(v(j).Name, 'Octave'), k = true; break, end
end

end
