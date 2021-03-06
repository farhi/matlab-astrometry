classdef astrometry < handle
  % astrometry: A Matlab class to annotate astrophotography images (identify objects/astrometry)
  %
  % #matlab-astrometry
  %
  % Purpose
  % -------
  %
  %  This Matlab class allows to use the astrometry.net software, either installed 
  %  locally, or through internet connection, in order to solve (annotate) 
  %  astrophotography images.
  %
  % Syntax/Usage
  % ------------
  %
  %  First navigate to the matlab-astrometry directory or type:
  %
  %  addpath /path/to/matlab-astrometry
  %
  %  Then use:
  %
  %  as = astrometry;
  %    Create a solver, but does not solve. 
  %    Use solve(as, file) or local(as, file) or web(as, file) afterwards.
  %
  %  as = astrometry(file, ...); image(as);
  %    Solve the given astrophotography image with local or web method. 
  %    Then plot the result. Additional arguments may include name/value pairs
  %    (see example below):
  %
  %      ra:      approximate RA coordinate  (e.g. deg or  'hh:mm:ss')
  %      dec:     approximate DEC coordinate (e.g. deg or 'deg:mm:ss')
  %      radius:  approximate field size     (in deg)
  %      scale-low:   lower estimate of the field coverage (in [deg], e.g. 0.1)
  %      scale-high:  upper estimate of the field coverage (in [deg], e.g. 180)
  %      object:  name of an object on field (string, e.g. 'M 16')
  %
  %  These two syntaxes will try first any local astrometry.net installation, and
  %  if failed, the http://nova.astrometry.net/ service.
  %
  %  When the annotation has ended, an 'annotationEnd' event is triggered. You
  %  may monitor this with e.g.:
  %    as = astrometry('examples/M13-2018-05-19.jpg');
  %    addlistener(as, 'annotationEnd', @(src,evt)disp('annotation just end'))
  %
  % Going further
  % -------------
  %
  %  as = LOAD(as, dir); IMAGE(as);
  %    Read an existing Astrometry.net set of files stored in a given directory.
  %    The directory may contain WCS, CORR, RDLS, JSON, and image.
  %    Then plot the result. This allows to get previous data files, or obtained
  %    externally, and label them. The 'as' astrometry object must have been used
  %    to solve or import astrometry data.
  %
  %  [x,y] = SKY2XY(as, ra, dec)
  %    Convert a RA/DEC set of coordinates (in [deg] or 'hh:mm:ss'/'deg::mm:ss')
  %    into pixel coordinates on the image. The 'as' astrometry object must have 
  %    been used to solve or import astrometry data.
  %
  %  [ra, dec] = XY2SKY(as, x,y)
  %  [ra, dec] = XY2SKY(as, x,y, 'string')
  %    Convert pixel coordinates on the image into a RA/DEC set of coordinates 
  %    (in [deg]). When given a 'string' argument, the result is given in 
  %    'hh:mm:ss'/'deg:mm:ss'. The 'as' astrometry object must have been used
  %    to solve or import astrometry data.
  %
  %  f = FINDOBJ(as,'object name')
  %    Return information about a named object (star, deep sky object) from the 
  %    data base. Example: astrometry.findobj('M33')
  %
  %  LOCAL(as, file, ...);
  %    Explicitly use the local 'solve-field' astrometry.net installation.
  %    See above for the additional arguments.
  %
  %  WEB(as, file, ...);
  %    Explicitly use the http://nova.astrometry.net/ web service.
  %    See above for the additional arguments.
  %
  %  WEB(as)
  %    For a solved image, the corresponding sky view is displayed on 
  %    http://www.sky-map.org . The 'as' astrometry object must have been used
  %    to solve or import astrometry data.
  %
  % Using results
  % -------------
  %   Once an image has been solved with the 'as' object, you can use the
  %   astrometry results.
  %
  %   The annotation is done asynchronously, and the Matlab prompt is recovered.
  %   You may use getstatus(as) to inquire for the solve-plate status
  %   (running, success, failed).
  %   To wait for the end of the annotation, use waitfor(as). stop(as) aborts it.
  %
  % - as.result.RA and as.result.Dec provide the center coordinates of the 
  %   field (in [deg]), while as.result.RA_hms and as.result.Dec_dms provide the 
  %   'HH:MM:SS' and 'Deg:MM:SS' coordinates. 
  % - The field rotation wrt sky is stored in as.result.rotation. 
  % - The pixel scale is given in [arcmin/pixel] as as.result.pixel_scale. 
  % - The field extension is given with its bounds as as.result.RA_min, as.result.RA_max,
  %   as.result.Dec_min, and as.result.Dec_min. 
  % - The constellation name is stored in as.result.Constellation.
  %
  % Improving the plate-solve efficiency
  % ------------------------------------
  %
  %  To facilitate the plate-solve/annotation of images, you may:
  %
  %  - specify the field size with additional arguments such as: 
  %     astrometry(..., 'scale-low', 0.5, 'scale-high',2)
  %
  %  - provide an initial guess for the location, and its range, such as:
  %     astrometry('examples/M13-2018-05-19.jpg', ...
  %       'ra','16:33:51','dec','30:39:35','radius', 2)
  %
  %  - provide the name of an object on field, such as:
  %     astrometry('examples/M13-2018-05-19.jpg','object','m 13','radius',2)
  % 
  %
  %  - add more star data bases (e.g. 2MASS over Tycho2).
  %
  % Examples
  % --------
  %
  %  as=astrometry('examples/M13-2018-05-19.jpg','scale-low', 0.5, 'scale-high',2);
  %  image(as);
  %
  % Methods
  % -------
  % - findobj   find a given object in catalogs.    
  % - getstatus return the astrometry status (success, failed)
  % - image     show the solve-plate image with annotations
  % - load      load astrometry files (WCS,FITS) from a directory   
  % - local     loads an image and identifies its objects using local solve-field  
  % - plot      show the solve-plate image with annotations. Same as image.   
  % - sky2xy    convert RA,Dec coordinates to x,y pixels on image   
  % - solve     solve an image field. Plot further results with IMAGE method.   
  % - stop      ends any current annotation and reset the object.   
  % - visible   return/display all visible objects on image   
  % - waitfor   waits for completion of the annotation   
  % - web       loads an image and identifies its objects using web service   
  % - xy2sky    convert pixel image coordinates to RA,Dec 
  %
  % Installation:
  % -------------
  %
  %  Local installation (recommended)
  %
  %    On Linux systems, install the 'astrometry.net' package, as well as the 
  %    'tycho2' data base. On Debian-class systems, this is achieved with:
  %
  %       sudo apt install astrometry.net astrometry-data-tycho2 sextractor
  %    
  %    On other systems, you will most probably need to compile it.
  %    See: http://astrometry.net/doc/build.html
  %    RedHat/Arch/MacOSX have specific installation instructions.
  %
  %    If you have images spanning on very tiny areas (e.g. much smaller than a 
  %    degree), you will most probably need to install the '2MASS' data base.
  %
  %  Using the web service.
  %
  %    You will need Python to be installed, and to have a 'NOVA astrometry API' key.
  %    Enter the API_KEY when prompt, or set it with:
  %
  %    as = astrometry;
  %    as.api_key = 'blah-blah';
  %    as.web(file, ...)
  %
  % Credit: 
  %
  %    sky2xy and xy2sky from E. Ofek http://weizmann.ac.il/home/eofek/matlab/
  %
  % (c) E. Farhi, 2018. GPL2.

  properties
  
    api_key    = '';  % api-key for nova.astrometry.net
     % example: 'kvfubnepntofzpcl' 'ghqpqhztzychczjh'
     % from: https://git.kpi.fei.tuke.sk/TP/ExplorationOfInterstellarObjects/blob/master/src/sk/tuke/fei/kpi/tp/eoio/AstrometryAPI.java

    result     = [];  % results from the annotation or empty when failed
    filename   = '';  % the image to annotate
    status     = 'init';  % can be: running, failed, success
    catalogs   = [];  % catalogs of common objects
    vargin     = [];  % arguments stored at instantiation for reuse
    autoplot   = false; % when true, display annotated image on success
    figure     = [];  % last figure shown
    duration   = 0;
  end % properties
  
  properties (Access=private)
    process_java = [];
    process_dir  = [];
    timer        = [];
    starttime    = [];
  end % private properties
  
  properties (Constant=true)
    executables  = find_executables;  % search for executables
  end % shared properties
  
  events
    annotationStart
    annotationEnd
    idle
    busy
  end
  
  methods
  
    function self=astrometry(filename, varargin)
      % ASTROMETRY Loads an image and identifies its objects using astrometry.net
      % 
      % as = astrometry;
      %   Create a solver, but does not solve. 
      %   Use local(as, file) or web(as, file) afterwards
      % as = astrometry(file, ...);
      %   Solve the given astrophotography image with local or web method.
      %
      % input(optional):
      %    filename: an image to annotate
      %    any name/value pair as:
      %      ra:      approximate RA coordinate  (e.g. deg or  'hh:mm:ss')
      %      dec:     approximate DEC coordinate (e.g. deg or 'deg:mm:ss')
      %      radius:  approximate field size     (in deg)
      %      scale-low:   lower estimate of the field coverage (in [deg], e.g. 0.1)
      %      scale-high:  upper estimate of the field coverage (in [deg], e.g. 180)
      %      object:  name of an object on field (string, e.g. 'M 16')
      % 
      % Example:
      %   as=astrometry('M33.jpg','scale-low', 0.5, 'scale-high',2);
      
      % handle input arguments
      removeme = [];
      for index=1:2:numel(varargin)
        if ischar(varargin{index})
          switch varargin{index}
          case 'catalogs'
            self.catalogs = varargin{index+1};
            removeme = [ removeme index index+1 ];
          case 'autoplot'
            self.autoplot = true;
            removeme = [ removeme index ];
          end
        end
      end
      if ~isempty(removeme)
        varargin(removeme) = [];
      end
      if ~isempty(varargin)
        self.vargin = varargin;
      end
      if isempty(self.catalogs)
        self.catalogs = getcatalogs;
      end
      
      if nargin
        % first try with the local plate solver
        [self.result, filename]      = self.solve(filename, 'solve-field', varargin{:});
        % if fails or not installed, use the web service
        if isempty(self.result) && ~isempty(filename)
          self.solve(filename, 'web', varargin{:});
        end
        % image(self);
      end
      
    end % astrometry
    
    function self = local(self, filename, varargin)
      % LOCAL Loads an image and identifies its objects using local solve-field
      %
      % as = local(astrometry, file, ...);
      %   Solve the given astrophotography image with local method.
      %
      % input(optional):
      %    filename: an image to annotate
      %    any name/value pair as:
      %      ra:      approximate RA coordinate  (e.g. deg or  'hh:mm:ss')
      %      dec:     approximate DEC coordinate (e.g. deg or 'deg:mm:ss')
      %      radius:  approximate field size     (in deg)
      %      scale-low:   lower estimate of the field coverage (in [deg], e.g. 0.1)
      %      scale-high:  upper estimate of the field coverage (in [deg], e.g. 180)
      %      object:  name of an object on field (string, e.g. 'M 16')
      % 
      % Example:
      %   as=local(astrometry, 'M33.jpg','scale-low', 0.5, 'scale-high',2);
      if isempty(varargin) && ~isempty(self.vargin), varargin = self.vargin; end
      [ret, filename] = solve(self, filename, 'solve-field', varargin{:});
    end % local (annotation)
    
    function self = web(self, filename, varargin)
      % WEB Loads an image and identifies its objects using web service
      %
      % as = web(astrometry, file, ...);
      %   Solve the given astrophotography image with web method.
      % web(as)
      %   Once solved, the field is displayed on http://www.sky-map.org
      %
      % input(optional):
      %    filename: an image to annotate
      %    any name/value pair as:
      %      ra:      approximate RA coordinate  (e.g. deg or  'hh:mm:ss')
      %      dec:     approximate DEC coordinate (e.g. deg or 'deg:mm:ss')
      %      radius:  approximate field size     (in deg)
      %      scale-low:   lower estimate of the field coverage (in [deg], e.g. 0.1)
      %      scale-high:  upper estimate of the field coverage (in [deg], e.g. 180)
      %      object:  name of an object on field (string, e.g. 'M 16')
      % 
      % Example:
      %   as=web(astrometry, 'M33.jpg','scale-low', 0.5, 'scale-high',2);
      if nargin == 1 && ~isempty(self.result) && isstruct(self.result)
        % display a Sky-Map.org view of the astrometry field
        sz = max([ self.result.RA_max-self.result.RA_min self.result.Dec_max-self.result.Dec_min ]);
        z  = 160.*2.^(0:-1:-8); % zoom levels in deg in sky-map
        z  = find(sz*4 > z, 1);
        if isempty(z), z=9; end
        url = sprintf([ 'http://www.sky-map.org/?ra=%f&de=%f&zoom=%d' ...
          '&show_grid=1&show_constellation_lines=1' ...
          '&show_constellation_boundaries=1' ...
          '&show_const_names=0&show_galaxies=1&img_source=DSS2' ], ...
          self.result.RA/15, self.result.Dec, z);
          % open in system browser
          open_system_browser(url);
      else
        if isempty(varargin) && ~isempty(self.vargin), varargin = self.vargin; end
        if nargin < 2 || isempty(filename)
          filename = self.filename;
        end
        [ret, filename] = solve(self, filename, 'web', varargin{:});
      end
    end % web
    
    function [ret, filename] = solve(self, filename, method, varargin)
      % SOLVE Solve an image field. Plot further results with IMAGE method.
      %
      %  solve(astrometry, filename, method, ...)
      %
      %  input:
      %    filename: an image to annotate
      %    method:   'solve-field' (default) or 'nova'
      %    any name/value pair as:
      %      ra:      approximate RA coordinate  (e.g. deg or  'hh:mm:ss')
      %      dec:     approximate DEC coordinate (e.g. deg or 'deg:mm:ss')
      %      radius:  approximate field size     (in deg)
      %      scale-low:   lower estimate of the field coverage (in [deg], e.g. 0.1)
      %      scale-high:  upper estimate of the field coverage (in [deg], e.g. 180)
      %      object:  name of an object on field (string, e.g. 'M 16')
      %      
      % example: as.solve('M33.jpg')
      %          as.solve('M33.jpg','default','ra','01:33:51','dec','30:39:35','radius', 2)
      %          as.solve('M33.jpg','default','ra','01:33:51','dec','30:39:35','radius', 2,...
      %             'scale-low',0.5,'scale-high',2)
      
      % varargin: ra, dec, radius, scale-low, scale-high
      
      ret = [];
      
      if nargin < 3,      method = ''; end
      if isempty(method), method = 'solve-field'; end
      
      if ishold(self)
        disp([ '[' datestr(now) ']: ' mfilename ': current solver is RUNNING.' ]);
        return; 
      end
      
      isnova = strcmp(method, 'nova')           || strcmp(method, 'web') ...
            || strcmp(method, 'astrometry.net') || strcmp(method, 'client.py');
            
      % check if executables are available      
      if isnova
        if isempty(self.executables.client_py) ...
        || isempty(self.executables.python), return; end
      else
        if isempty(self.executables.solve_field), return; end
      end
      
      % request image if missing
      if nargin < 2, filename = ''; end
      if isempty(filename), filename = self.filename; end
      if iscell(filename)
        ret={};
        for index=1:numel(filename)
          ret{end+1} = solve(self, filename{index}, method, varargin{:});
        end
        return
      end
      if ischar(filename) && isempty(filename)
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
        if isequal(filename,0), filename = ''; return; end
        filename = fullfile(pathname, filename);
      elseif isempty(filename)
        return; % nothing to do
      end
      % when given as RGB
      if isnumeric(filename)
        rgb = filename;
        filename = tempname;
        imwrite(im, [ filename '.png' ], 'png', ...
          'Author', getenv('USER'), ...
          'CreationTime', datestr(now), ...
          'Software', mfilename);
      end
      % handle distant images
      if ischar(filename) && any(strncmp(filename, {'http','ftp:'}, 4))
        [p,f,e] = fileparts(filename);
        tmpfile = [ tempname e ];
        filename = urlwrite(filename, tmpfile);
      end
      if ~ischar(filename), filename = ''; return;
      elseif any(strcmpi(filename,{'/dev/null','null','none'}))
        filename = []; 
        return
      elseif isempty(dir(filename))
        disp([ '[' datestr(now) ']: ' mfilename ': ERROR: invalid file name "' filename '"' ]);
        filename = []; 
        return
      end
      self.filename = filename;
      
      % build the command line
      if isempty(self.process_dir) % first call
        self.process_dir = tempname;
      elseif ~isempty(dir(self.process_dir))
        % clean any previous annotation
        rmdir(self.process_dir, 's');
      end
      if ~isdir(self.process_dir)
        mkdir(self.process_dir);
      end
      d=self.process_dir;

      if isnova
        % is there an API_KEY ? request it if missing...
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
      
        cmd = [ self.executables.python ' ' self.executables.client_py ' --wait' ];
        if ~isempty(self.api_key)
          cmd = [ cmd ' --apikey=' self.api_key ];
        end
        cmd = [ cmd ' --upload='   filename ];
        cmd = [ cmd ' --annotate=' fullfile(d, 'results.json') ];
        cmd = [ cmd ' --newfits='  fullfile(d, 'results.fits') ];
        cmd = [ cmd ' --kmz='      fullfile(d, 'results.kml') ];
        cmd = [ cmd ' --corr='     fullfile(d, 'results.corr') ];
      else
        cmd = [ self.executables.solve_field  ];
        cmd = [ cmd  ' '           filename ];
        if ~isempty(self.executables.sextractor)
          % highly improves annotation efficiency
          cmd = [ cmd ' --use-sextractor' ];
        end
        cmd = [ cmd ' --dir '      d ];
        cmd = [ cmd ' --new-fits ' fullfile(d, 'results.fits') ];
        cmd = [ cmd ' --rdls '     fullfile(d, 'results.rdls') ];
        cmd = [ cmd ' --corr '     fullfile(d, 'results.corr') ' --tag-all' ];
        if ~isempty(self.executables.wcs2kml)
          cmd = [ cmd ' --kmz '      fullfile(d, 'results.kml') ' --no-tweak' ];
        end
      end
      cmd = [ cmd ' --wcs='      fullfile(d, 'results.wcs') ];

      % handle arguments: target name/findobj
      if nargin > 3
        obj = []; rmarg = [];
        for f=1:numel(varargin)
          if ischar(varargin{f}) && f < numel(varargin) ...
            && any(strcmp(varargin{f}, {'findobj','goto','obj','object'})) ...
            if ischar(varargin{f+1})
              obj = self.findobj(varargin{f+1});
            elseif isstruct(varargin{f+1})
              obj = varargin{f+1};
            end
            rmarg = [ f f+1 ];
          elseif isstruct(varargin{f})
            obj = varargin{f};
            rmarg = f;
          end
        end
        if ~isempty(obj) && isstruct(obj) ...
          && isfield(obj, 'RA') && isfield(obj, 'DEC')
          varargin(rmarg) = [];
          varargin{end+1} = 'ra';           varargin{end+1} = obj.RA;
          varargin{end+1} = 'dec';          varargin{end+1} = obj.DEC;
          if isfield(obj,'NAME')
            disp([ '[' datestr(now) ']: ' mfilename ': goto "' obj.NAME '"' ]);
          end
        end
      end
      
      % handle additional arguments in name/value pairs
      if nargin > 3 && mod(numel(varargin), 2) == 0
        for f=1:2:numel(varargin)
          if ischar(varargin{f})
            if isnova
              if     strcmp(varargin{f}, 'scale-low')  varargin{f}='scale-lower'; 
              elseif strcmp(varargin{f}, 'scale-high') varargin{f}='scale-upper'; end
            else
              if     strcmp(varargin{f}, 'scale-lower') varargin{f}='scale-low'; 
              elseif strcmp(varargin{f}, 'scale-upper') varargin{f}='scale-high'; end
            end
            cmd = [ cmd ' --' varargin{f} '=' num2str(varargin{f+1}) ];
          end
        end
      end

      % execute command
      disp(cmd)
      self.status   = 'running';
      self.starttime= clock;
      notify(self, 'busy');
      disp([ '[' datestr(now) ']: ' mfilename ': ' method ' please wait (may take e.g. few minutes)...' ])
      
      % create the timer for auto update
      if isempty(self.timer) || ~isa(self,timer, 'timer') || ~isvalid(self.timer)
        self.timer  = timer('TimerFcn', @TimerCallback, ...
          'Period', 5.0, 'ExecutionMode', 'fixedDelay', 'UserData', self, ...
          'Name', mfilename);
      end
      if strcmp(self.timer.Running, 'off') start(self.timer); end
          
      % launch a Java asynchronous command
      self.process_java = java.lang.Runtime.getRuntime().exec(cmd);
      
      % we shall monitor the completion with a timer
      notify(self, 'annotationStart');
      ret = cmd;
    end % solve
    
    function ret = load(self, d, varargin)
      % LOAD Load astrometry files (WCS,FITS) from a directory
      %   LOAD(astrometry, directory);
      %   The directory may contain WCS, CORR, RDLS or JSON, and image.
      %   No solve plate is performed, only data is read.
      %
      %   LOAD(astrometry, image, ...) starts solve-plate annotation, same as SOLVE.
      if nargin < 2, d = self.process_dir; end
      if ~isdir(d)
        % loading an image ?
        try
          im = imread(d);
          solve(self, d, varargin{:});
          return
        end
      end
      
      self.result = getresult(d, self);
      if isempty(self.result)
        self.status = 'failed';
        self.process_java = [];
      else
        self.status = 'success';
        % is the image available ? use one from the directory
        if ~exist(self.filename, 'file')
          % search in the result directory
          d = [ dir(fullfile(self.result.dir, '*.png')) ; 
                dir(fullfile(self.result.dir, '*.fits')) ];
          if ~isempty(d), 
            d=d(1);
            self.filename = fullfile(self.result.dir, d.name);
          end
        end
      end
      ret = self.result;
      
    end % load
    
    function ret=getstatus(self)
      % GETSTATUS Return the astrometry status (success, failed)
      ret=self.status;
    end % getstatus
    
    function st = ishold(self)
      % ISHOLD get the solver status (IDLE, BUSY)
      %   st = ISHOLD(s) returns 1 when the solver is BUSY.
      st = ~isempty(self.process_java);
    end % ishold
    
    function fig = plot(self, varargin)
      % PLOT Show the solve-plate image with annotations. Same as image.
      %
      %   as=astrometry(file);
      %   plot(as);
      fig = image(self, varargin{:});
    end % plot
    
    function [fig] = image(self, mag)
      % IMAGE Show the solve-plate image with annotations
      %   IMAGE(as, mag) limits the objects up to given magnitude.
      %
      %   as=astrometry(file);
      %   image(as);

      if nargin < 2, mag = inf; end
      
      fig = [];
      if ~ischar(self.filename) || isempty(self.filename) || isempty(dir(self.filename)), return; end
      try
        im  = imread(self.filename);
      catch ME
        getReport(ME)
        disp([ '[' datestr(now) ']: ' mfilename  ': ERROR: can not read image ' self.filename ]);
        return
      end 
      fig = figure('Name', [ mfilename ': ' self.filename ]);
      self.figure = fig;
      image(im);
      clear im;
      
      ret = self.result;
      
      % set title
      [p,f,e] = fileparts(self.filename);
      if isfield(ret, 'Constellation')
        title([ f e ' in ' ret.Constellation ]);
      else
        title([ f e ]);
      end
      
      % overlay results
      if ~isempty(self.result) && strcmp(self.status, 'success') && isfield(self.result, 'RA_hms')
        
        hold on
        % central coordinates
        sz = self.result.size/2;
        h  = plot(sz(1), sz(2), 'r+'); set(h, 'MarkerSize', 16);
        hcmenu = uicontextmenu;
        uimenu(hcmenu, 'Label', '<html><b>Field center</b></html>');
        uimenu(hcmenu, 'Label', [ 'RA=  ' ret.RA_hms ]);
        uimenu(hcmenu, 'Label', [ 'DEC= ' ret.Dec_dms ]);
        uimenu(hcmenu, 'Label', [ 'Rotation= ' num2str(ret.rotation) ' [deg]' ]);
        set(h, 'UIContextMenu', hcmenu);
        
        % get list of visible objects
        v = visible(self, mag);
        
        % make sure we do not plot too many
        if numel(v) > 1500
          v = v(1:1500);
        end
        
        for index=1:numel(v)
          % find all objects from data base within bounds
          this = v(index);

          % stars in green, DSO in cyan
          if strcmp(this.catalog,'stars'), c = 'g'; sz = 12;
          else                             c = 'c'; 
          end
          x = this.X; y = this.Y;
          
          % plot symbol
          if isfinite(this.SIZE) && this.SIZE > 5
            h = plot(x,y, [ c 'o' ]); 
            set(h, 'MarkerSize', ceil(this.SIZE));
          else
            h = plot(x,y, [ c 's' ]); this.SIZE=12;
          end
          
          % context menu
          hcmenu = uicontextmenu;            
          uimenu(hcmenu, 'Label', [ 'RA=  ' this.RA ' [' num2str(this.RA_deg) ' deg]' ]);
          uimenu(hcmenu, 'Label', [ 'DEC= ' this.DEC ' [' num2str(this.DEC_deg) ' deg]' ]);
          uimenu(hcmenu, 'Label', [ '<html><b>' this.NAME '</html></b>' ], 'Separator','on');
          uimenu(hcmenu, 'Label', [ 'TYPE: ' this.TYPE ]);
          if isfinite(this.MAG) && this.MAG > 0
            uimenu(hcmenu, 'Label', [ 'MAGNITUDE= ' num2str(this.MAG)  ]);
          end
          if isfinite(this.DIST) && this.DIST > 0
            uimenu(hcmenu, 'Label', [ 'DIST= ' sprintf('%.3g', this.DIST*3.262) ' [ly]' ]);
          end
          set(h, 'UIContextMenu', hcmenu);
          if numel(v) < 1500
            t=text(x+this.SIZE,y-this.SIZE,this.NAME); set(t,'Color', c);
          end
        end % for

      end % success
      
    end % image
    
    function [x,y] = sky2xy(self, ra,dec)
      % SKY2XY Convert RA,Dec coordinates to x,y pixels on image
      %
      % input:
      %   ra,dec: RA and Dec [deg]
      % output:
      %   x,y:    pixel coordinates
      x = []; y = [];
      if isempty(self.result), return; end
      if ~isscalar(ra)
        ra = getra(ra);
      end
      if ~isscalar(dec)
        dec = getdec(dec);
      end
      [x,y] = sky2xy_tan(self.result.wcs.meta, ...
         ra*pi/180, dec*pi/180);                         % MAAT Ofek (private)
    end
    
    function [ra,dec] = xy2sky(self, x,y, str)
      % XY2SKY Convert pixel image coordinates to RA,Dec
      %
      % input:
      %   x,y:    pixel coordinates
      %   str: when 
      % output:
      %   ra,dec: RA and Dec [deg]
      ra = []; dec = [];
      if isempty(self.result), return; end
      
      if nargin > 3, str=true; else str=false; end
      [ra, dec] = xy2sky_tan(self.result.wcs.meta, x,y); % MAAT Ofek (private)
      ra =ra *180/pi;
      dec=dec*180/pi;
      if str
        ra = getra( ra/15,  'string');
        dec= getdec(dec, 'string');
      end
    end
    
    function found = findobj(self, name)
      % FINDOBJ Find a given object in catalogs. 
      catalogs = fieldnames(self.catalogs);
      found = [];
      
      % check first for name without separator
      if ~any(name == ' ')
        [n1,n2]  = strtok(name, '0123456789');
        found = findobj(self, [ n1 ' ' n2 ]);
        if ~isempty(found) return; end
      end
      namel= strtrim(lower(name));
      for f=catalogs(:)'
        catalog = self.catalogs.(f{1});
        if ~isfield(catalog, 'MAG'), continue; end
        NAME = lower(catalog.NAME);
        NAME = regexprep(NAME, '\s*',' ');
        % search for name
        index = find(~cellfun(@isempty, strfind(NAME, [ ';' namel ';' ])));
        if isempty(index)
        index = find(~cellfun(@isempty, strfind(NAME, [ namel ';' ])));
        end
        if isempty(index)
        index = find(~cellfun(@isempty, strfind(NAME, [ ';' namel ])));
        end
        if isempty(index)
        index = find(~cellfun(@isempty, strfind(NAME, [ namel ])));
        end
        if ~isempty(index)
          found.index   = index(1);
          found.catalog = f{1};
          found.RA      = catalog.RA(found.index);
          found.DEC     = catalog.DEC(found.index);
          found.MAG     = catalog.MAG(found.index);
          found.TYPE    = catalog.TYPE{found.index};
          found.NAME    = catalog.NAME{found.index};
          found.DIST    = catalog.DIST(found.index);
          break;
        end
      end

      if ~isempty(found)
        disp([ '[' datestr(now) ']: ' mfilename  ': Found object ' name ' as: ' found.NAME ])
        if found.DIST > 0
          disp(sprintf('  %s: Magnitude: %.1f ; Type: %s ; Dist: %.3g [ly]', ...
            found.catalog, found.MAG, found.TYPE, found.DIST*3.262 ));
        else
          disp(sprintf('  %s: Magnitude: %.1f ; Type: %s', ...
            found.catalog, found.MAG, found.TYPE ));
        end
        if ~isempty(self.result) && isfield(self.result, 'RA_hms')
          ret = self.result;
          if  (ret.RA_min <= found.RA ...
            &               found.RA  <= ret.RA_max ...
            & ret.Dec_min<= found.DEC ...
            &               found.DEC <= ret.Dec_max)
            disp([ '    ' found.NAME ' is within the image field.' ]);
          end
        end
      else
        disp([ '[' datestr(now) ']: ' mfilename  ': object ' name ' was not found.' ])
      end
    end % findobj
    
    function v = visible(self, mag)
      % VISIBLE Return/display all visible objects on image
      %   VISIBLE(as, mag) limits the objects up to given magnitude.
      
      v = [];
      if nargin < 2, mag = inf; end

      if ~isempty(self.result) && isfield(self.result, 'RA_hms')
        
        ret = self.result;
        
        for catalogs = {'stars','deep_sky_objects'}
          % find all objects from data base within bounds
          catalog = self.catalogs.(catalogs{1});
          
          found = find(...
              ret.RA_min <= catalog.RA ...
            &               catalog.RA  <= ret.RA_max ...
            & ret.Dec_min<= catalog.DEC ...
            &               catalog.DEC <= ret.Dec_max ...
            & catalog.MAG <= mag);
          
          if ~isempty(found)
            disp([ '[' datestr(now) ']: ' mfilename  ': ' num2str(numel(found)) ...
              ' referenced ' catalogs{1} ' in field.' ])
          end
          for index=1:numel(found)
            obj     = found(index);
            ra      = catalog.RA(obj);
            dec     = catalog.DEC(obj);
            [x,y]   = self.sky2xy(ra, dec);
            
            % ignore when not on image
            if x < 1 || x > ret.size(1) || y < 1 || y > ret.size(2), continue; end
            
            this.RA = getra(ra/15, 'string');
            this.DEC= getdec(dec,  'string');
            this.NAME=catalog.NAME{obj};
            this.TYPE=catalog.TYPE{obj};
            this.MAG =catalog.MAG(obj);
            this.SIZE=catalog.SIZE(obj); % arcmin
            this.DIST=catalog.DIST(obj);
            this.X   =x;
            this.Y   =y;
            this.catalog=catalogs{1};
            this.RA_deg = ra;
            this.DEC_deg= dec;
            
            if isempty(v), v = this;
            else v(end+1) = this; end
            
          end % object
          
        end % catalogs
        
        % sort by increasing magnitude (decreasing brightness)
        [~,index]=sort([ v.MAG ]);
        v = v(index);
        
        % display the list
        if nargout == 0
          disp(self.filename)
          disp 'TYPE            MAG  RA              DEC                 DIST  NAME'
          disp '----------------------------------------------------------------------'
          for index=1:numel(v)
            this = v(index);
            % display the list
            fprintf(1, '%-12s  %5.1f  %-14s  %-14s  %8.2g  %s\n', ...
              this.TYPE, this.MAG, this.RA, this.DEC, this.DIST*3.262, this.NAME);
          end
        end
        
      end % success
      
    end % visible
    
    function disp(self)
      % DISP Display Astrometry object (details)
      
      if ~isempty(inputname(1))
        iname = inputname(1);
      else
        iname = 'ans';
      end
      if isdeployed || ~usejava('jvm') || ~usejava('desktop'), id=class(self);
      else id=[  '<a href="matlab:doc ' class(self) '">' class(self) '</a> ' ...
                 '(<a href="matlab:methods ' class(self) '">methods</a>,' ...
                 '<a href="matlab:image(' iname ');">plot</a>,' ...
                 '<a href="matlab:visible(' iname ');">visible...</a>)' ];
      end
      fprintf(1,'%s = %s for "%s":\n',iname, id, self.filename);
      if ~isempty(self.result) && strcmp(self.status, 'success') && isfield(self.result, 'RA_hms')
        % get list of visible objects, extract brightest star and DSO
        v = visible(self); 
        mag = [ inf, inf ];  % star, dso
        name= { '',  ''};
        
        for index=1:numel(v)
          % find all objects from data base within bounds
          this = v(index);
          % stars in green, DSO in cyan
          if strcmp(this.catalog,'stars') obj_index=1;
          else                            obj_index=2;
          end
          if this.MAG < mag(obj_index)
            mag(obj_index) = this.MAG;
            name{obj_index}= [ this.NAME ' ' this.TYPE ];
          end
        end
        if isfield(self.result, 'Constellation')
          disp([ '  Constellation: ' self.result.Constellation ]);
        end
        disp([   '  RA:            ' self.result.RA_hms  ' [h:min:s]; ' ...
          num2str(self.result.RA)  ' [deg]']);
        disp([   '  DEC:           ' self.result.Dec_dms ' [deg:min:s]; ' ...
          num2str(self.result.Dec) ' [deg]' ]);
        for index=1:2
          if ~isempty(name{index})
            disp([ '    magnitude ' num2str(mag(index)) ' ' name{index} ])
          end
        end
        disp([   '  Rotation:      ' num2str(self.result.rotation) ' [deg] (to get sky view)' ]);
        disp([   '  Pixel scale:   ' num2str(self.result.pixel_scale) ' [arcsec/pixel]' ]);
        if isdeployed || ~usejava('jvm') || ~usejava('desktop')
          disp([ '  Results are in ' self.process_dir ]);
        else
          disp([ '  Results are in <a href="' self.process_dir '">' self.process_dir '</a>' ]);
        end
        builtin('disp',self)
        disp([ iname '.result:' ])
        disp(self.result);
      elseif isdir(self.process_dir)
        if isdeployed || ~usejava('jvm') || ~usejava('desktop')
          disp([ '  ' upper(self.status) ' in ' self.process_dir ]);
        else
          disp([ '  ' upper(self.status) ' in <a href="' self.process_dir '">' self.process_dir '</a>' ]);
        end
      else
        disp([ '  ' upper(self.status) ': use annotate(as,''filename'').' ]);
      end
    
    end % disp
    
    function display(self)
      % DISPLAY Display Astrometry object (short)
      
      if ~isempty(inputname(1))
        iname = inputname(1);
      else
        iname = 'ans';
      end
      if isdeployed || ~usejava('jvm') || ~usejava('desktop'), id=class(self);
      else id=[  '<a href="matlab:doc ' class(self) '">' class(self) '</a> ' ...
                 '(<a href="matlab:methods ' class(self) '">methods</a>,' ...
                 '<a href="matlab:image(' iname ');">plot</a>,' ...
                 '<a href="matlab:disp(' iname ');">more...</a>)' ];
      end
      if ishold(self)
        fprintf(1,'%s = %s for "%s" BUSY\n',iname, id, self.filename);
      else fprintf(1,'%s = %s for "%s"\n',iname, id, self.filename);
      end
    end % display
    
    function waitfor(self)
      % WAITFOR Waits for completion of the annotation
      while ishold(self)
        pause(5);
      end
    end % waitfor
    
    function stop(self)
    % STOP Ends any current annotation and reset the object.
      % clear the timer
      if ~isempty(self.timer) && isa(self.timer,'timer')
        stop(self.timer); delete(self.timer); 
      end
      self.timer=[];
      p = self.process_java;
      if ~isempty(p) && isjava(p)
        try
          p.destroy;
          disp([ '[' datestr(now) ']: ' mfilename  ': abort current annotation...' ])
        end
      end
      self.process_java = [];
      self.status = 'failed';
    end % stop
    
    function c = get_catalogs(self)
      % GET_CATALOGS Get the loaded catalogs
      c = self.catalogs;
    end
    
    function st = get_state(self)
      % GET_STATE Return the astrometry state, e.g. BUSY, FAILED, SUCCESS.
      st = self.status;
    end
    
  end % methods
  
end % astrometry class

% ------------------------------------------------------------------------------

function ret = getresult(d, self)
  % getresult: extract WCS and star matching information from the output files.
  %
  % input:
  %   d: directory where astrometry.net results are stored.
  
  ret = [];
  for file={'results.wcs','wcs.fits'}
    if exist(fullfile(d, file{1}), 'file')
      ret.wcs  = read_fits(fullfile(d, file{1}));
      
      % get image center and print it
      if isfield(ret.wcs,'meta') && isfield(ret.wcs.meta,'CRVAL1')
        wcs = ret.wcs.meta;
        ret.wcs.meta.CD = [ wcs.CD1_1 wcs.CD1_2 ; 
                            wcs.CD2_1 wcs.CD2_2 ];
                            
        % get central coordinates
        
        ret.size     = [ wcs.IMAGEW wcs.IMAGEH ];
        sz  = ret.size/2;

        [ret.RA, ret.Dec] = xy2sky_tan(ret.wcs.meta, sz(1), sz(2)); % MAAT Ofek (private)
        ret.RA=ret.RA*180/pi; ret.Dec=ret.Dec*180/pi;
        ret.RA_hms   = getra(ret.RA/15,   'string');
        ret.Dec_dms  = getdec(ret.Dec, 'string');
        % compute pixel scale
        ret.pixel_scale = sqrt(abs(wcs.CD1_1 * wcs.CD2_2  - wcs.CD1_2 * wcs.CD2_1))*3600; % in arcsec/pixel
        % compute rotation angle
        ret.rotation = atan2(wcs.CD2_1, wcs.CD1_1)*180/pi;
                            
        % compute RA,Dec image bounds
        RA=[]; Dec=[];
        [RA(end+1), Dec(end+1)]     = xy2sky_tan(ret.wcs.meta, wcs.IMAGEW, wcs.IMAGEH);
        [RA(end+1), Dec(end+1)]     = xy2sky_tan(ret.wcs.meta, 1         , wcs.IMAGEH);
        [RA(end+1), Dec(end+1)]     = xy2sky_tan(ret.wcs.meta, wcs.IMAGEW, 1);
        [RA(end+1), Dec(end+1)]     = xy2sky_tan(ret.wcs.meta, 1         , 1);
        RA = RA*180/pi; Dec = Dec*180/pi;
        ret.RA_min  = min(RA);
        ret.RA_max  = max(RA);
        ret.Dec_min = min(Dec);
        ret.Dec_max = max(Dec);
        
        % identify constellation we are in
        [m, index] = min( (ret.RA - self.catalogs.constellations.RA).^2 ...
                        + (ret.Dec- self.catalogs.constellations.DEC).^2 );
        ret.Constellation = self.catalogs.constellations.Name{index};
      end
    end
  end % for
  for file={'results.rdls','rdls.fits'}
    if exist(fullfile(d, file{1}), 'file')
      ret.rdls = read_fits(fullfile(d, file{1}));
    end
  end
  for file={'results.corr','corr.fits'}
    if exist(fullfile(d, file{1}), 'file')
      ret.corr = read_fits(fullfile(d, file{1}));
    end
  end
  if exist(fullfile(d, 'results.json'), 'file')
    ret.json = loadjson(fullfile(d, 'results.json'));
  end
  if ~isempty(ret)
    ret.dir     = d;
  end
  
end % getresult

% ------------------------------------------------------------------------------
function executables = find_executables
  % find_executables: locate executables, return a structure
  
  % stored here so that they are not searched for further calls
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
  for exe =  { 'solve-field','sextractor','python','python3', 'wcs2kml', 'client.py' }
  
    for try_target={ ...
      fullfile(this_path, [ exe{1} ext ]), ...
      fullfile(this_path, [ exe{1} ]), ...
      [ exe{1} ext ], ... 
      exe{1} }
      
      if exist(try_target{1}, 'file')
        status = 0; result = 'OK';
      else
        [status, result] = system([ precmd try_target{1} ' --version' ]); % run from Matlab
      end
      
      name = strrep(exe{1}, '-','_');
      name = strrep(name  , '.','_');

      if status ~= 127
        % the executable is found
        if strcmp(name, 'python3') name = 'python'; end
        executables.(name) = try_target{1};
        disp([ '[' datestr(now) ']: ' mfilename ': found ' exe{1} ' as ' try_target{1} ])
        break
      else
        executables.(name) = [];
      end
    end
  
  end
  
  found_executables = executables;
  
end % find_executables

function [ra_h, ra_min, ra_s] = getra(ra, str)
  % getra: convert any input RA (in hours) into h and min
  
  if nargin > 1, str = true; else str=false; end
  if ischar(ra)
    ra = repradec(ra);
  end
  if isstruct(ra) && isfield(ra, 'RA')
    ra = ra.RA;
  end
  if isstruct(ra)
    if isfield(ra, 'h'),   ra_h   = ra.h; end
    if isfield(ra, 'min'), ra_min = ra.min; end
    if isfield(ra, 's'),   ra_s   = ra.s; end
  end
  if isnumeric(ra)
    if isscalar(ra)
      ra_h   = fix(ra); 
      ra_min = abs(ra - ra_h)*60; 
      ra_s   = abs(ra_min - fix(ra_min))*60;
      ra_min = fix(ra_min);
    elseif numel(ra) == 2
      ra_h   = ra(1);     
      ra_min = abs(ra(2)); 
      ra_s   = abs(ra_min - fix(ra_min))*60;
      ra_min = fix(ra_min);
    elseif numel(ra) == 3
      ra_h   = ra(1);
      ra_min = ra(2);
      ra_s   = ra(3);
    end
  else
    disp([ '[' datestr(now) ']: ' mfilename  ': invalid RA.' ])
    disp(ra)
  end
  if nargout == 1
    if str
      ra_h = [ num2str(ra_h)    ':' num2str(ra_min)  ':' num2str(ra_s) ];
    else
      ra_h = ra_h+ra_min/60 + ra_s/3600;
    end
  elseif nargout == 2
    ra_min = ra_min + ra_s/60;
  end
end % getra

function [dec_deg, dec_min, dec_s] = getdec(dec, str)
  % getdec: convert any input DEC into deg and min
  
  if nargin > 1, str = true; else str=false; end
  if ischar(dec)
    dec = repradec(dec);
  end
  if isstruct(dec) && isfield(dec, 'DEC')
    dec = dec.DEC;
  end
  if isstruct(dec)
    if isfield(dec, 'deg') dec_deg = dec.deg; end
    if isfield(dec, 'min') dec_min = dec.min; end
    if isfield(dec, 'min') dec_s   = dec.s; end
  end
  if isnumeric(dec)
    if isscalar(dec)
      dec_deg = floor(dec); 
      dec_min = abs(dec - dec_deg)*60;
      dec_s   = abs(dec_min - fix(dec_min))*60;
      dec_min = fix(dec_min);
    elseif numel(dec) == 2
      dec_deg = dec(1);   
      dec_min = abs(dec(2));
      dec_s   = abs(dec_min - fix(dec_min))*60;
      dec_min = fix(dec_min);
    elseif numel(dec) == 3
      dec_deg = dec(1);   
      dec_min = dec(2);
      dec_s   = dec(3);
    end
  else
    disp([ '[' datestr(now) ']: ' mfilename ': invalid DEC' ])
    disp(dec)
  end
  if nargout == 1
    if str
      dec_deg = [ num2str(dec_deg) ':' num2str(dec_min) ':' num2str(dec_s) ];
    else
      dec_deg = dec_deg + dec_min/60 + dec_s/3600;
    end
  elseif nargout == 2
    dec_min = dec_min + dec_s/60;
  end
end % getdec

function str = repradec(str)
  % repradec: replace string stuff and get it into num
  str = lower(str);
  for rep = {'h','m','s',':','°','deg','d','''','"'}
    str = strrep(str, rep{1}, ' ');
  end
  str = str2num(str);
end % repradec

% ------------------------------------------------------------------------------

function catalogs = getcatalogs
  % getcatalogs: load catalogs for stars and DSO.

  % stored here so that they are not loaded for further calls
  persistent loaded_catalogs  
  
  if ~isempty(loaded_catalogs)
    catalogs = loaded_catalogs; 
    return
  end
  
  % load catalogs: objects, stars
  disp([ '[' datestr(now) ']: ' mfilename  ': Welcome ! Loading Catalogs:' ]);
  catalogs = load(mfilename);
  
  % display available catalogs
  for f=fieldnames(catalogs)'
    name = f{1};
    if ~isempty(catalogs.(name))
      num  = numel(catalogs.(name).RA);
      if isfield(catalogs.(name), 'Description')
        desc = catalogs.(name).Description;
      else desc = ''; end
      disp([ mfilename ': ' name ' with ' num2str(num) ' entries.' ]);
      disp([ '  ' desc ])
    end
  end

  loaded_catalogs = catalogs;
end % getcatalogs

% ------------------------------------------------------------------------------

function ret=open_system_browser(url)
  % opens URL with system browser. Returns non zero in case of error.
  if strncmp(url, 'file://', length('file://'))
    url = url(8:end);
  end
  ret = 1;
  if ismac,      precmd = 'DYLD_LIBRARY_PATH= ;';
  elseif isunix, precmd = 'LD_LIBRARY_PATH= ; '; 
  else           precmd=''; end
  if ispc
    ret=system([ precmd 'start "' url '"']);
  elseif ismac
    ret=system([ precmd 'open "' url '"']);
  else
    [ret, message]=system([ precmd 'xdg-open "' url '"']);
  end
end % open_system_browser

% ------------------------------------------------------------------------------

function TimerCallback(src, evnt)
  % TimerCallback: update status/view from timer event
  self = get(src, 'UserData');
  if isvalid(self)
    try
      % check if any astrometry job is running
      exitValue = 0;
      if ~isempty(self.process_java) && isjava(self.process_java)
        try
          exitValue = self.process_java.exitValue; % will raise error if process still runs
          active  = 0;
        catch ME
          % still running
          if isempty(self.process_java) || ~isjava(self.process_java)
            active  = 0;
          else
            active  = 1;
          end
        end
        % not active anymore: process has ended.
        if ~active
          if exitValue ~= 0
            self.result = []; 
            self.status = 'failed';
          else
            load(self);
            if ~isempty(self.result) self.status = 'success';
            else self.status = 'failed'; end
          end
          disp([ '[' datestr(now) ']: ' mfilename ': annotation end: ' upper(self.status) ' for ' self.filename '. exit value=' num2str(exitValue) ]);
          % clear the timer
          stop(src); delete(src); self.timer=[];
          self.process_java = [];
          self.duration = etime(clock, self.starttime);
          
          notify(self, 'annotationEnd');
          notify(self, 'idle');
          if self.autoplot
            image(self);
            assignin('base','ans', self);
            ans = self;
            disp(ans);
          end
          beep;
        end
      end
    catch ME
      getReport(ME)
    end
    
  else delete(src); end
  
end % TimerCallback

