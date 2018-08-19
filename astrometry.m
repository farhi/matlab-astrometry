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
      % cmd = [ cmd ' --newfits='  fullfile(d, 'results.fits') ];
      % cmd = [ cmd ' --kmz='      fullfile(d, 'results.kml') ];

      disp(cmd)
      disp([ mfilename ': please wait (may take e.g. few minutes)...' ])
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
      end
    
    end % client
    
    function [ret, filename] = solve(self, filename, varargin)
      %  input:
      %    filename: an image to annotate
      %    any name/value pairs as:
      %      ra:      approximate RA coordinate  (e.g. deg or 'hh:mm:ss')
      %      dec:     approximate DEC coordinate (e.g. deg or 'deg:mm:ss')
      %      radius:  approximate field size      (in deg)
      %      kmz:     output KML file             (requires wcs2kml)
      %      depth:   number of field objects to look at (e.g. 1-50)
      
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

      % --dir tempdir
      % solve-field IC5146-Cocon-2018-08-11.jpg  --ra 21:53:24 --dec 47:16:14 --radius 2  -l 600 -L 0.5 -H 2 --overwrite
    
    end % solve
    
    
  end % methods
  
end % astrometry

% ------------------------------------------------------------------------------
function executables = find_executables
  % find_executables: locate executables, return a structure
  
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
  
end % find_executables
