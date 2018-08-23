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

meta = [];
% search for some 'Keywords'
for extname=fieldnames(s)'
  if isfield(s.(extname{1}), 'Keywords')
    kw = s.(extname{1}).Keywords;
    for index=1:size(kw,1)
      name = kw{index,1};
      val  = kw{index,2};
      if isempty(val)
        val = kw{index,3};
      end
      if ischar(val)
        val = strtrim(val); 
      end
      if ~isfield(meta, name)
        meta.(name) = val;
      else
        if ~iscell(meta.(name)), meta.(name) = { meta.(name) };
        else meta.(name) = [ meta.(name) ; val ]; end
      end
    end
  end
end

s.meta = meta;

% get data from the BinaryTable.
% The fields are described in the BinaryTable.Keywords
% and the data  itself  is in the data.bintable as a cell with arrays
% the TTYPE* fields indicate what this is, and the order in the table

if isfield(meta, 'TFIELDS') && isfield(meta, 'XTENSION') ...
  && strcmp(meta.XTENSION, 'BINTABLE') && numel(s.data.bintable) == meta.TFIELDS
  f = fieldnames(meta);
  for index=1:meta.TFIELDS
    name = [ 'TTYPE' num2str(index) ];
    if isfield(meta, name)
      name = meta.(name);
      val  = s.data.bintable{index};
      s.data.(name) = val;
    end
  end  
  
end
