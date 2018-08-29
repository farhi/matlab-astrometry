classdef astrometry < handle
  % astrometry: A Matlab class to annotate astrophotography images (identify objects/astrometry)
  %
  %  Methods:
  %  
  %   as=astrometry(filename)
  %   image(as)
  %   load(astrometry, dir)
  %   annotation(astrometry, filename, ...)
  %   web(astrometry, filename, ...)
  %   sky2sx(as, ra, dec)
  %   xy2sky(as, x, y)
  %
  % Credit: sky2xy and xy2sky from E. Ofek http://weizmann.ac.il/home/eofek/matlab/
  % (c) E. Farhi, 2018. GPL2.

  properties
  
    api_key    = '';  % api-key for nova.astrometry.net
     % example: 'kvfubnepntofzpcl' 'ghqpqhztzychczjh'
     % from: https://git.kpi.fei.tuke.sk/TP/ExplorationOfInterstellarObjects/blob/master/src/sk/tuke/fei/kpi/tp/eoio/AstrometryAPI.java
    executables= [];  % executables
    result     = [];
    filename   = [];
    status     = 'init';  % can be: running, failed, success
    log        = '';
    catalogs   = [];
    
  end % properties
  
  methods
  
    function self=astrometry(filename, varargin)
      % astrometry: loads an image and identifies its objects using astrometry.net
      % 
      % as = astrometry;
      %   Create a solver, but does not solve. 
      %   Use annotate(as, file) or web(as, file) afterwards
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
      % 
      % Example:
      %   as=astrometry('M33.jpg','scale-low', 0.5, 'scale-high',2);
      self.executables = find_executables;
      self.catalogs    = getcatalogs;
      
      if nargin
        % first try with the local plate solver
        [self.result, filename]      = self.solve(filename, 'solve-field', varargin{:});
        % if fails or not installed, use the web service
        if isempty(self.result)
          self.solve(filename, 'web', varargin{:});
        end
        % image(self);
      end
      
    end % astrometry
    
    function self = annotation(self, filename, varargin)
      % astrometry.annotation: loads an image and identifies its objects using local solve-field
      %
      % as = annotation(astrometry, file, ...);
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
      % 
      % Example:
      %   as=annotation(astrometry, 'M33.jpg','scale-low', 0.5, 'scale-high',2);
      [ret, filename] = solve(self, filename, 'solve-field', varargin{:});
    end
    
    function self = web(self, filename, varargin)
      % astrometry.web: loads an image and identifies its objects using web service
      %
      % as = web(astrometry, file, ...);
      %   Solve the given astrophotography image with web method.
      %
      % input(optional):
      %    filename: an image to annotate
      %    any name/value pair as:
      %      ra:      approximate RA coordinate  (e.g. deg or  'hh:mm:ss')
      %      dec:     approximate DEC coordinate (e.g. deg or 'deg:mm:ss')
      %      radius:  approximate field size     (in deg)
      %      scale-low:   lower estimate of the field coverage (in [deg], e.g. 0.1)
      %      scale-high:  upper estimate of the field coverage (in [deg], e.g. 180)
      % 
      % Example:
      %   as=web(astrometry, 'M33.jpg','scale-low', 0.5, 'scale-high',2);
      [ret, filename] = solve(self, filename, 'web', varargin{:});
    end
    
    function [ret, filename] = solve(self, filename, method, varargin)
      % astrometry.solve: solve an image field. Plot further results with image method.
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
      %      
      % example: as.solve('M33.jpg')
      %          as.solve('M33.jpg','default','ra','01:33:51','dec','30:39:35','radius', 2)
      %          as.solve('M33.jpg','default','ra','01:33:51','dec','30:39:35','radius', 2,...
      %             'scale-low',0.5,'scale-high',2)
      
      % varargin: ra, dec, radius, scale-low, scale-high
      
      ret = [];
      
      if nargin < 3,      method = ''; end
      if isempty(method), method = 'solve-field'; end
      
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
      self.filename = filename;
      
      % build the command line
      if ismac,      precmd = 'DYLD_LIBRARY_PATH= ;';
      elseif isunix, precmd = 'LD_LIBRARY_PATH= ; '; 
      else           precmd=''; end
      
      d = tempname;
      if ~isdir(d), mkdir(d); end

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
      
        cmd = [ precmd self.executables.python ' ' self.executables.client_py ' --wait' ];
        if ~isempty(self.api_key)
          cmd = [ cmd ' --apikey=' self.api_key ];
        end
        cmd = [ cmd ' --upload='   filename ];
        cmd = [ cmd ' --annotate=' fullfile(d, 'results.json') ];
        cmd = [ cmd ' --newfits='  fullfile(d, 'results.fits') ];
        cmd = [ cmd ' --kmz='      fullfile(d, 'results.kml') ];
      else
        cmd = [ precmd self.executables.solve_field  ];
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
      
      % handle additional arguments in name/value pairs
      if nargin > 3 && mod(numel(varargin), 2) == 0
        for f=1:2:numel(varargin)
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

      % execute command
      disp(cmd)
      t0 = clock;
      self.status   = 'running';
      disp([ mfilename ': ' method ': please wait (may take e.g. few minutes)...' ])
      [status, self.log] = system(cmd);
      
      if status ~= 0
        self.result = []; 
        self.status = 'failed';
      else
        load(self, d);
      end
      
      if ~isempty(self.result)
        ret.duration= etime(clock, t0);
      end
      disp([ mfilename ': ' upper(self.status) ': plate solve for image ' filename ' using ' method ])
    
    end % solve
    
    function load(self, d)
      % astrometry.load: load astrometry files (WCS,FITS) from a directory
      %
      %   The directory may contain WCS, CORR, RDLS or JSON, and image.
      %   No solve plate is performed, only data is read.
      %
      % as = load(astrometry, directory);
      self.result = getresult(d, self);
      if isempty(self.result)
        self.status = 'failed';
      else
        self.status = 'success';
        % is the image available ? use one from the directory
        if ~exist(self.filename, 'file')
          % search in the result directory
          d = dir(fullfile(self.result.dir, '*.png'));
          if ~isempty(d), 
            d=d(1);
            self.filename = fullfile(self.result.dir, d.name);
          end
        end
      end
      
    end
    
    function fig = image(self)
      % astrometry.image: show the solve-plate image with annotations
      %
      %   as=astrometry(file);
      %   image(as);
      
      fig = [];
      try
        im  = imread(self.filename);
      catch ME
        getReport(ME)
        disp([ mfilename ': ERROR: can not read image ' self.filename ]);
        return
      end 
      fig = figure('Name', [ mfilename ': ' self.filename ]);
      image(im); title(self.filename);
      im_sz = size(im);
      clear im;
      
      % overlay results
      if ~isempty(self.result) && strcmp(self.status, 'success')
        ret = self.result;
        hold on
        
        % identify constellation we are in
        
        
        % central coordinates
        sz = im_sz/2; % [ height width layers ]
        h  = plot(sz(2), sz(1), 'r+'); set(h, 'MarkerSize', 16);
        hcmenu = uicontextmenu;
        uimenu(hcmenu, 'Label', '<html><b>Field center</b></html>');
        uimenu(hcmenu, 'Label', [ 'RA=  ' ret.RA_hms ]);
        uimenu(hcmenu, 'Label', [ 'DEC= ' ret.Dec_dms ]);
        set(h, 'UIContextMenu', hcmenu);
        
        for catalogs = {'stars','deep_sky_objects'}
          % find all objects from data base within bounds
          catalog = self.catalogs.(catalogs{1});
          
          found = find(...
              ret.RA_min <= catalog.RA ...
            &               catalog.RA  <= ret.RA_max ...
            & ret.Dec_min<= catalog.DEC ...
            &               catalog.DEC <= ret.Dec_max);
          
          for index=1:numel(found)
            obj     = found(index);
            ra      = catalog.RA(obj);
            dec     = catalog.DEC(obj);
            [x,y]   = self.sky2xy(ra, dec);
            
            % ignore when not on image
            if x < 1 || x > im_sz(2) || y < 1 || y > im_sz(1), continue; end
            
            ra_hms  = getra(ra/15, 'string');
            dec_dms = getdec(dec,  'string');
            name    = catalog.NAME{obj};
            typ     = catalog.TYPE{obj};
            mag     = catalog.MAG(obj);
            sz      = catalog.SIZE(obj); % arcmin
            dist    = catalog.DIST(obj);

            % stars in green, DSO in cyan
            if strcmp(catalogs{1},'stars'), c = 'g'; sz = 12;
            else                            c = 'c'; 
            end
            
            % plot symbol
            if isfinite(sz) && sz > 5
              h = plot(x,y, [ c 'o' ]); 
              set(h, 'MarkerSize', ceil(sz));
            else
              h = plot(x,y, [ c 's' ]); sz=12;
            end
            
            % context menu
            hcmenu = uicontextmenu;            
            uimenu(hcmenu, 'Label', [ 'RA=  ' ra_hms ]);
            uimenu(hcmenu, 'Label', [ 'DEC= ' dec_dms ]);
            uimenu(hcmenu, 'Label', [ '<html><b>' name '</html></b>' ], 'Separator','on');
            uimenu(hcmenu, 'Label', [ 'TYPE: ' typ ]);
            if isfinite(mag) && mag > 0
              uimenu(hcmenu, 'Label', [ 'MAGNITUDE= ' num2str(mag)  ]);
            end
            if isfinite(dist) && dist > 0
              uimenu(hcmenu, 'Label', [ 'DIST= ' sprintf('%.3g', dist) ' [ly]' ]);
            end
            set(h, 'UIContextMenu', hcmenu);
            t=text(x+sz,y-sz,name); set(t,'Color', c);
          end % object
        end % catalogs

      end % success
      
    end % image
    
    function [x,y] = sky2xy(self, ra,dec)
      % astrometry.sky2xy: convert RA,Dec coordinates to x,y pixels on image
      %
      % input:
      %   ra,dec: RA and Dec [deg]
      % output:
      %   x,y:    pixel coordinates
      x = []; y = [];
      if isempty(self.result), return; end
      [x,y] = sky2xy_tan(self.result.wcs.meta, ...
         ra*pi/180, dec*pi/180);                         % MAAT Ofek (private)
    end
    
    function [ra,dec] = xy2sky(self, x,y)
      % astrometry.xy2sky: convert pixel image coordinates to RA,Dec
      %
      % input:
      %   x,y:    pixel coordinates
      % output:
      %   ra,dec: RA and Dec [deg]
      ra = []; dec = [];
      if isempty(self.result), return; end
      [ra, dec] = xy2sky_tan(self.result.wcs.meta, x,y); % MAAT Ofek (private)
      ra =ra *180/pi;
      dec=dec*180/pi;
    end
    
  end % methods
  
end % astrometry

% ------------------------------------------------------------------------------

function ret = getresult(d, self)
  % getresult: extract WCS and star matching information from the output files.
  %
  % input:
  %   d: directory where astrometry.net results are stored.
  
  if isdeployed || ~usejava('jvm') || ~usejava('desktop')
    disp([ '  Results are in: ' d ]);
  else
    disp([ '  Results are in: <a href="' d '">' d '</a>' ]);
  end
  
  ret = [];
  if exist(fullfile(d, 'results.wcs'), 'file')
    ret.wcs  = read_fits(fullfile(d, 'results.wcs'));
    
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
      
      disp([ 'Field center:    RA =' ret.RA_hms  ' [h:min:s]    ; ' ...
        num2str(ret.RA)  ' [deg]']);
      disp([ '                 DEC=' ret.Dec_dms ' [deg:min:s]  ; ' ...
        num2str(ret.Dec) ' [deg]' ]);
      % compute pixel scale
      ret.pixel_scale = sqrt(abs(wcs.CD1_1 * wcs.CD2_2  - wcs.CD1_2 * wcs.CD2_1))*3600; % in arcsec/pixel
      disp([ 'Pixel scale:         ' num2str(ret.pixel_scale) ' [arcsec/pixel]' ]);
      % compute rotation angle
      ret.rotation = atan2(wcs.CD2_1, wcs.CD1_1)*180/pi;
      disp([ 'Field rotation:      ' num2str(ret.rotation) ' [deg] (to get sky view)' ]);
                          
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
    end
  end
  if exist(fullfile(d, 'results.rdls'), 'file')
    ret.rdls = read_fits(fullfile(d, 'results.rdls'));
  end
  if exist(fullfile(d, 'results.corr'), 'file')
    ret.corr = read_fits(fullfile(d, 'results.corr'));
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
    disp([ mfilename ': invalid RA.' ])
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
    disp([ mfilename ': invalid DEC' ])
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
  disp([ mfilename ': Welcome ! Loading Catalogs:' ]);
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

