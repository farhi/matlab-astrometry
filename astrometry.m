classdef astrometry < handle
  % astrometry: A Matlab class to annotate astrophotography images (identify objects/astrometry)
  %
  %  Purpose
  %
  % (c) E. Farhi, 2018. GPL2.

  properties
  
    api_key    = '';  % api-key for nova.astrometry.net
     % example: 'kvfubnepntofzpcl' 'ghqpqhztzychczjh'
     % from: https://git.kpi.fei.tuke.sk/TP/ExplorationOfInterstellarObjects/blob/master/src/sk/tuke/fei/kpi/tp/eoio/AstrometryAPI.java
    executables= [];  % executables
    
  end % properties
  
  methods
  
    function self=astrometry(filename, varargin)
      % astrometry: loads an image and identifies its objects using astrometry.net
      % 
      self.executables = find_executables;
      
    end % astrometry
    
    function [ret, filename] = client(self, filename, varargin)
      %  input:
      %    filename: an image to annotate
      
      % varargin: ra, dec, radius, scale-lower, scale-upper
      
      ret = [];
      if isempty(self.executables.client_py) ...
      || isempty(self.executables.python), return; end
      
      if nargin < 2, filename = ''; end
      if isempty(filename)
        % request an image file to solve
        [filename, pathname] = uigetfile( ...
          {'*.JPG;*.JPEG;*.jpg;*.jpeg;*.GIF;*.gif;*.PNG;*.png;*.FITS;*.fits;*.FTS;*.fts', ...
              'All supported images (JPG,PNG,GIF,FITS)';
           '*.JPG;*.JPEG;*.jpg;*.jpeg', 'JPEG image (JPG)';
           '*.GIF;*.gif', 'GIF image (GIF)';
           '*.PNG;*.png', 'PNG image (PNG)';
           '*.FITS;*.fits;*.FTS;*.fts','FITS image (FITS)';
           '*.*',  'All Files (*.*)'}, ...
           [ mfilename ': Pick an astrophotography image to solve' ]);
        if isequal(filename,0), return; end
        filename = fullfile(pathname, filename);
      end
      
      if isempty(self.api_key) && isempty(getenv('AN_API_KEY'))
        % request the API_KEY via a dialogue
        op.Resize      = 'on';
        op.WindowStyle = 'normal';   
        op.Interpreter = 'tex';
        NL = sprintf('\n');
        prompt = [ ...
          '{\color{blue}Enter a nova.astrometry.net API key}' NL ...
          '  (e.g. "slsratwckfdyxhjq")' NL ...
          'Create an account at {\color{red}http://nova.astrometry.net} (free)' NL ...
          'Connect to your account and click on the "{\color{blue}API}" tab to get the key' NL ...
          'or define the {\color{blue}AN\_API\_KEY} environment variable' ];
        answer = inputdlg( prompt, 'Input Astrometry.net API key', 1, { 'slsratwckfdyxhjq' }, op);
        if isempty(answer), return; end
        self.api_key = answer{1};
      end
      
      if ismac,      precmd = 'DYLD_LIBRARY_PATH= ;';
      elseif isunix, precmd = 'LD_LIBRARY_PATH= ; '; 
      else           precmd=''; end
      
      d = tempname;
      if ~isdir(d), mkdir(d); end
      
      cmd = [ precmd self.executables.python ' ' self.executables.client_py ' --wait' ];
      if ~isempty(self.api_key)
        cmd = [ cmd ' --apikey=' self.api_key ];
      end
      cmd = [ cmd ' --upload='   filename ];
      cmd = [ cmd ' --wcs='      fullfile(d, 'results.wcs') ];
      cmd = [ cmd ' --annotate=' fullfile(d, 'results.json') ];
      cmd = [ cmd ' --newfits='  fullfile(d, 'results.fits') ];
      cmd = [ cmd ' --kmz='      fullfile(d, 'results.kml') ];
      
      % handle additional arguments in name/value pairs
      if nargin > 3 && mod(numel(varargin), 2) == 0
        for f=1:2:numel(varargin)
          cmd = [ cmd ' --' varargin{f} '=' num2str(varargin{f+1}) ];
        end
      end

      disp(cmd)
      t0 = clock;
      disp([ mfilename ': client.py: please wait (may take e.g. few minutes)...' ])
      [status, result] = system(cmd);
      
      if status ~= 0
        disp(result)
        error([ mfilename ': FAILED: plate solve for image ' filename ' using http://nova.astrometry.net' ])
      end
      disp([ mfilename ': SUCCESS: plate solve for image ' filename ' using http://nova.astrometry.net' ])
      disp([ '  Results are in: ' d ]);
      
      if exist(fullfile(d, 'results.json'), 'file')
        ret.json = loadjson(fullfile(d, 'results.json'));
      end
      if exist(fullfile(d, 'results.wcs'), 'file')
        ret.wcs  = read_fits(fullfile(d, 'results.wcs'));
        % get image center and print it
        if isfield(ret.wcs,'meta') && isfield(ret.wcs.meta,'CRVAL1')
          ret.ra  = ret.wcs.meta.CRVAL1;
          [ret.ra_h, ret.ra_min, ret.ra_s] = getra(ret.ra/15);
          ret.dec  = ret.wcs.meta.CRVAL2;
          [ret.dec_deg, ret.dec_min, ret.dec_s] = getdec(ret.dec);
          disp([ '  RA=  ' num2str(ret.ra_h)    ':' num2str(ret.ra_min)  ':' num2str(ret.ra_s) ])
          disp([ '  DEC= ' num2str(ret.dec_deg) ':' num2str(ret.dec_min) ':' num2str(ret.dec_s) ])
        end
      end
      ret.dir     = d;
      ret.duration= etime(clock, t0);
      
    
    end % client
    
    function [ret, filename] = solve(self, filename, varargin)
      %  input:
      %    filename: an image to annotate
      %    any name/value pair as:
      %      ra:      approximate RA coordinate  (e.g. deg or 'hh:mm:ss')
      %      dec:     approximate DEC coordinate (e.g. deg or 'deg:mm:ss')
      %      radius:  approximate field size      (in deg)
      %      scale-low:   lower estimate of the field coverage (in [deg], e.g. 0.1)
      %      scale-high:  upper estimate of the field coverage (in [deg], e.g. 180)
      %      kmz:     output KML file             (requires wcs2kml)
      %      depth:   number of field objects to look at (e.g. 1-50)
      %      
      % example: as.solve('M33.jpg')
      %          as.solve('M33.jpg','ra','01:33:51','dec','30:39:35','radius', 2)
      %          as.solve('M33.jpg','ra','01:33:51','dec','30:39:35','radius', 2,...
      %             'scale-low',0.5,'scale-high',2)
      
      % varargin: ra, dec, radius, scale-low, scale-high
      
      ret = [];
      if isempty(self.executables.solve_field), return; end
      
      if nargin < 2, filename = ''; end
      if isempty(filename)
        % request an image file to solve
        [filename, pathname] = uigetfile( ...
          {'*.JPG;*.JPEG;*.jpg;*.jpeg;*.GIF;*.gif;*.PNG;*.png;*.FITS;*.fits;*.FTS;*.fts', ...
              'All supported images (JPG,PNG,GIF,FITS)';
           '*.JPG;*.JPEG;*.jpg;*.jpeg', 'JPEG image (JPG)';
           '*.GIF;*.gif', 'GIF image (GIF)';
           '*.PNG;*.png', 'PNG image (PNG)';
           '*.FITS;*.fits;*.FTS;*.fts','FITS image (FITS)';
           '*.*',  'All Files (*.*)'}, ...
           [ mfilename ': Pick an astrophotography image to solve' ]);
        if isequal(filename,0), return; end
        filename = fullfile(pathname, filename);
      end
      
      if ismac,      precmd = 'DYLD_LIBRARY_PATH= ;';
      elseif isunix, precmd = 'LD_LIBRARY_PATH= ; '; 
      else           precmd=''; end
      
      d = tempname;
      if ~isdir(d), mkdir(d); end
      
      cmd = [ precmd self.executables.solve_field  ];
      cmd = [ cmd  ' '           filename ];
      if ~isempty(self.executables.sextractor)
        % highly improves annotation efficiency
        cmd = [ cmd ' --use-sextractor' ];
      end
      cmd = [ cmd ' --dir '      d ];
      cmd = [ cmd ' --wcs '      fullfile(d, 'results.wcs') ];
      cmd = [ cmd ' --new-fits ' fullfile(d, 'results.fits') ];
      cmd = [ cmd ' --rdls '     fullfile(d, 'results.rdls') ];
      cmd = [ cmd ' --corr '     fullfile(d, 'results.corr') ' --tag-all' ];
      if ~isempty(self.executables.wcs2kml)
        cmd = [ cmd ' --kmz '      fullfile(d, 'results.kml') ' --no-tweak' ];
      end
      
      % handle additional arguments in name/value pairs
      if nargin > 3 && mod(numel(varargin), 2) == 0
        for f=1:2:numel(varargin)
          cmd = [ cmd ' --' varargin{f} '=' num2str(varargin{f+1}) ];
        end
      end

      disp(cmd)
      t0 = clock;
      disp([ mfilename ': solve-field: please wait (may take e.g. few minutes)...' ])
      [status, result] = system(cmd);
      
      if status ~= 0
        disp(result)
        error([ mfilename ': FAILED: plate solve for image ' filename ' using solve-field' ])
      end
      disp([ mfilename ': SUCCESS: plate solve for image ' filename ' using solve-field' ])
      disp([ '  Results are in: ' d ]);
      
      if exist(fullfile(d, 'results.wcs'), 'file')
        ret.wcs  = read_fits(fullfile(d, 'results.wcs'));
        % get image center and print it
        if isfield(ret.wcs,'meta') && isfield(ret.wcs.meta,'CRVAL1')
          ret.ra  = ret.wcs.meta.CRVAL1;
          [ret.ra_h, ret.ra_min, ret.ra_s] = getra(ret.ra/15);
          ret.dec  = ret.wcs.meta.CRVAL2;
          [ret.dec_deg, ret.dec_min, ret.dec_s] = getdec(ret.dec);
          disp([ '  RA=  ' num2str(ret.ra_h)    ':' num2str(ret.ra_min)  ':' num2str(ret.ra_s) ])
          disp([ '  DEC= ' num2str(ret.dec_deg) ':' num2str(ret.dec_min) ':' num2str(ret.dec_s) ])
        end
      end
      if exist(fullfile(d, 'results.rdls'), 'file')
        ret.rdls = read_fits(fullfile(d, 'results.rdls'));
      end
      if exist(fullfile(d, 'results.corr'), 'file')
        ret.corr = read_fits(fullfile(d, 'results.corr'));
      end
      ret.dir     = d;
      ret.duration= etime(clock, t0);
    
    end % solve
    
    
  end % methods
  
end % astrometry

% ------------------------------------------------------------------------------
function executables = find_executables
  % find_executables: locate executables, return a structure
  
  persistent found_executables
  
  if ~isempty(found_executables)
    executables = found_executables;
    return
  end
  
  if ismac,      precmd = 'DYLD_LIBRARY_PATH= ;';
  elseif isunix, precmd = 'LD_LIBRARY_PATH= ; '; 
  else           precmd=''; end
  
  if ispc, ext='.exe'; else ext=''; end
  
  executables = [];
  this_path   = fullfile(fileparts(which(mfilename)));
  
  % what we may use
  for exe =  { 'solve-field','sextractor','python', 'wcs2kml', 'client.py' }
  
    for try_target={ [ exe{1} ext ], exe{1}, ...
      fullfile(this_path, [ exe{1} ]), ... 
      fullfile(this_path, [ exe{1} ext ])}
      
      if exist(try_target{1}, 'file')
        status = 0; result = 'OK';
      else
        [status, result] = system([ precmd try_target{1} ' --version' ]); % run from Matlab
      end
      
      name = strrep(exe{1}, '-','_');
      name = strrep(name  , '.','_');

      if status ~= 127
        % the executable is found
        executables.(name) = try_target{1};
        disp([ mfilename ': found ' exe{1} ' as ' try_target{1} ])
        break
      else
        executables.(name) = [];
      end
    end
  
  end
  
  found_executables = executables;
  
end % find_executables

function [ra_h, ra_min, ra_s] = getra(ra)
  % getra: convert any input RA into h and min
  ra_h = []; ra_min = [];
  if ischar(ra)
    ra = repradec(ra);
  end
  if isstruct(ra) && isfield(ra, 'RA')
    ra = ra.RA;
  elseif isstruct(ra) && isfield(ra, 'h') && isfield(ra, 'min')
    ra_h = ra.h;
    ra_min = ra.min;
    return
  end
  if isnumeric(ra)
    if isscalar(ra)
      ra_h   = fix(ra); ra_min = abs(ra - ra_h)*60;
    elseif numel(ra) == 2
      ra_h = ra(1);     ra_min = abs(ra(2));
    elseif numel(ra) == 3
      ra_h = ra(1);     ra_min = abs(ra(2))+abs(ra(3)/60);
    end
    if nargout == 3
      ra_s   = abs(ra_min - fix(ra_min))*60;
      ra_min = fix(ra_min);
    end
  else
    disp([ mfilename ': invalid RA.' ])
    disp(ra)
  end
  if nargout == 1
    ra_h = ra_h+ra_min/60;
  end
end % getra

function [dec_deg, dec_min, dec_s] = getdec(dec)
  % getdec: convert any input DEC into deg and min
  if ischar(dec)
    dec = repradec(dec);
  end
  if isstruct(dec) && isfield(dec, 'DEC')
    dec = dec.DEC;
  elseif isstruct(dec) && isfield(dec, 'deg') && isfield(dec, 'min')
    dec_deg = dec.deg;
    dec_min = dec.min;
    return
  end
  if isnumeric(dec)
    if isscalar(dec)
      dec_deg = floor(dec); dec_min = abs(dec - dec_deg)*60;
    elseif numel(dec) == 2
      dec_deg = dec(1);   dec_min = abs(dec(2));
    elseif numel(dec) == 3
      dec_deg = dec(1);   dec_min = abs(dec(2))+abs(dec(3)/60);
    end
    if nargout == 3
      dec_s   = abs(dec_min - fix(dec_min))*60;
      dec_min = fix(dec_min);
    end
  else
    disp([ mfilename ': invalid DEC' ])
    disp(dec)
  end
  if nargout == 1
    dec_deg = dec_deg + dec_min/60;
  end
end % getdec

function str = repradec(str)
  % repradec: replace string stuff and get it into num
  str = lower(str);
  for rep = {'h','m','s',':','°','deg','d','''','"'}
    str = strrep(str, rep{1}, ' ');
  end
  str = str2num(str);
end
