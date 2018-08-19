function s = read_fits(filename)
% read_fits Wrapper to fitsinfo/fitsread which reconstructs the FITS structure
%   s = read_fits(filename)
%
% (c) E.Farhi, ILL. License: EUPL.
% See also: read_edf, read_adsc, read_edf, read_sif, read_mar, read_spe, read_cbf, read_image, read_hbin

try
  s      = fitsinfo(filename);
catch
  s=[]; return;
end
% we read all FITS 'extnames'
s.data = [];
for extname={'primary','table','bintable','image','unknown'}
  try
    s.data.(extname{1}) = fitsread(filename, extname{1});
  end
end

