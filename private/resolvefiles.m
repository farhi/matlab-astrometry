function a = resolvefiles(a, t)
% resolvefiles: determine the fully qualified path of given names.
%   also resolves the wildcards expansion.
%
% output:
%   cellstr of filenames
  if nargin <2, t='Pick file(s)'; end
  if isempty(a)
    formats = imformats;
    
    % assemble the list of formats for the file selector
    f = {};
    f{end+1} = '*.*';
    f{end+1} = 'All Files (*.*)';
    f{end+1} = '*.jpg;*.JPG;*.jpeg;*.JPEG';
    f{end+1} = 'JPEG images (*.jpg)';
    f{end+1} = '*.png;*.PNG';
    f{end+1} = 'PNG images (*.jpg)';
    f{end+1} = '*.tif;*.TIF;*.tiff;*.TIFF';
    f{end+1} = 'TIFF images (*.tiff)';
    f{end+1} = '*.fit;*.FIT;*.fits;*.FITS;*.fts;*.FTS';
    f{end+1} = 'FITS images (*.fits)';
    if exist('readraw')
      f{end+1} = [...
       '*.raw;*.crw;*.cr2;*.kdc;*.dcr;*.mrw;*.arw;' ...
       '*.nef;*.nrw;*.dng;*.orf;*.ptx;*.pef;*.rw2;*.srw;*.raf;*.kdc' ...
       '*.RAW;*.CRW;*.CR2;*.KDC;*.DCR;*.MRW;*.ARW;' ...
       '*.NEF;*.NRW;*.DNG;*.ORF;*.PTX;*.PEF;*.RW2;*.SRW;*.RAF;*.KDC'];
      f{end+1} = 'RAW Camera Format (RAW)';
    end
    for index=1:numel(formats)
      this = formats(index);
      F=upper(this.ext);
      f{end+1} = [ sprintf('*.%s;', this.ext{:}) sprintf('*.%s;', F{:}) ];
      f{end+1} = this.description;
    end
    f = reshape(f, [ 2 numel(f)/2])';
   
    % get the files
    [a, pathname] = uigetfile( ...
          f, ...
          t, ...
          'MultiSelect', 'on');
    if isequal(a,0) || isempty(a), return; end
    a = strcat(pathname, a);
  else

    if ischar(a),    a = cellstr(a); end
    if iscell(a)
      new_a = {};
      for index=1:numel(a)
        this = a{index};
        if ischar(this) && strncmp(this,'file://',7)
          this = this(8:end);
        end
        if ischar(this)
          p    = fileparts(this);
          if isdir(this), p = this; end
          if isempty(p), p=pwd; end
          b    = dir(this); % in case we have a wildcard
          for n=1:numel(b)
            if ~b(n).isdir
              new_a{end+1} = fullfile(p, b(n).name);
            end
          end
        elseif isnumeric(this)
          new_a{end+1} = this;
        end
      end
      a = new_a;
    end
    
  end
end % resolvefiles
