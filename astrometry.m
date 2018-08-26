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
    result     = [];
    filename   = [];
    status     = '';  % can be: running, failed, success
    log        = '';
    catalogs   = [];
    
  end % properties
  
  methods
  
    function self=astrometry(filename, varargin)
      % astrometry: loads an image and identifies its objects using astrometry.net
      % 
      self.executables = find_executables;
      self.catalogs    = getcatalogs;
      
      if nargin
        % first try with the local plate solver
        [self.result, filename]      = self.solve(filename, varargin{:});
        % if fails or not installed, use the web service
        if isempty(self.result)
          self.client(filename, varargin{:});
        end
      end
      
    end % astrometry
    
    function [ret, filename] = client(self, filename, varargin)
      %  input:
      %    filename: an image to annotate
      %    any name/value pair as:
      %      ra:      approximate RA coordinate  (e.g. deg or 'hh:mm:ss')
      %      dec:     approximate DEC coordinate (e.g. deg or 'deg:mm:ss')
      %      radius:  approximate field size      (in deg)
      %      scale-lower: lower estimate of the field coverage (in [deg], e.g. 0.1)
      %      scale-upper: upper estimate of the field coverage (in [deg], e.g. 180)
      
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
        if     strcmp(varargin{f}, 'scale-low')  varargin{f}='scale-lower'; 
        elseif strcmp(varargin{f}, 'scale-high') varargin{f}='scale-upper';
        end
        for f=1:2:numel(varargin)
          cmd = [ cmd ' --' varargin{f} '=' num2str(varargin{f+1}) ];
        end
      end

      disp(cmd)
      t0 = clock;
      self.status   = 'running';
      self.filename = filename;
      disp([ mfilename ': client.py: please wait (may take e.g. few minutes)...' ])
      [status, self.log] = system(cmd);
      
      if status ~= 0
        disp([ mfilename ': FAILED: plate solve for image ' filename ' using http://nova.astrometry.net' ]);
        self.status = 'failed';
        return
      end
      disp([ mfilename ': SUCCESS: plate solve for image ' filename ' using http://nova.astrometry.net' ])
      
      ret = getresult(d);
      ret.duration= etime(clock, t0);
      self.status = 'success';
      
      self.result = ret;
    
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
          if     strcmp(varargin{f}, 'scale-lower') varargin{f}='scale-low'; 
          elseif strcmp(varargin{f}, 'scale-upper') varargin{f}='scale-high';
          end
          cmd = [ cmd ' --' varargin{f} '=' num2str(varargin{f+1}) ];
        end
      end

      disp(cmd)
      t0 = clock;
      self.status   = 'running';
      self.filename = filename;
      disp([ mfilename ': solve-field: please wait (may take e.g. few minutes)...' ])
      [status, self.log] = system(cmd);
      
      if status ~= 0
        disp([ mfilename ': FAILED: plate solve for image ' filename ' using solve-field' ])
        self.status = 'failed';
        return
      end
      disp([ mfilename ': SUCCESS: plate solve for image ' filename ' using solve-field' ])
      
      ret = getresult(d);
      ret.duration= etime(clock, t0);
      self.status = 'success';
      
      self.result = ret;
    
    end % solve
    
    function fig = image(self)
    
      fig = [];
      try
        im  = imread(self.filename);
      catch
        return
      end 
      fig = figure('Name', [ mfilename ': ' self.filename ]);
      image(im);
      sz = size(im);
      clear im;
      
      % overlay results
      if ~isempty(self.result) && strcmp(self.status, 'success')
        ret = self.result;
        hold on
        
        % central coordinates
        h = plot(sz(2)/2, sz(1)/2, 'y+'); set(h, 'MarkerSize', 16);
        hcmenu = uicontextmenu;
        uimenu(hcmenu, 'Label', '<html><b>Field centre</b></html>');
        uimenu(hcmenu, 'Label', [ 'RA= ' num2str(ret.ra_h)    ':' num2str(ret.ra_min)  ':' num2str(ret.ra_s) ]);
        uimenu(hcmenu, 'Label', [ 'DEC= ' num2str(ret.dec_deg) ':' num2str(ret.dec_min) ':' num2str(ret.dec_s) ]);
        set(h, 'UIContextMenu', hcmenu);
        % reference stars
        if isfield(ret, 'corr') && isfield(ret.corr.data, 'field_x')
          for index=1:numel(ret.corr.data.field_x)
            x   = ret.corr.data.field_x(index);
            y   = ret.corr.data.field_y(index);
            ra  = ret.corr.data.field_ra(index);
            dec = ret.corr.data.field_dec(index);
            [ra_h, ra_min, ra_s]      = getra(ra/15);
            [dec_deg, dec_min, dec_s] = getdec(dec);
            % identify the star names
            [minradec, st] = min( abs(ra - self.catalogs.stars.RA) + abs(dec - self.catalogs.stars.DEC) );
            
            if minradec < 1e-1
              name = self.catalogs.stars.NAME{st};
              typ  = self.catalogs.stars.TYPE{st};
              dist = self.catalogs.stars.DIST(st);
              mag  = self.catalogs.stars.MAG(st);
              h = plot(x,y,'go'); set(h, 'MarkerSize', 12);
              hcmenu = uicontextmenu;
              
              uimenu(hcmenu, 'Label', [ 'RA=  ' num2str(ra_h)    ':' num2str(ra_min)  ':' num2str(ra_s) ]);
              uimenu(hcmenu, 'Label', [ 'DEC= ' num2str(dec_deg) ':' num2str(dec_min) ':' num2str(dec_s) ]);
              uimenu(hcmenu, 'Label', [ '<html><b>' name '</html></b>' ], 'Separator','on');
              
              uimenu(hcmenu, 'Label', [ 'TYPE: ' typ ]);
              uimenu(hcmenu, 'Label', [ 'MAGNITUDE= ' num2str(mag)  ]);
              uimenu(hcmenu, 'Label', [ 'DIST= ' num2str(dist*3.262/1000, 2) ' [kly]' ]);
              set(h, 'UIContextMenu', hcmenu);
            end
            
          end
          
          % find all DSO within the reference star area
          min_ra = min(ret.corr.data.field_ra);
          max_ra = max(ret.corr.data.field_ra);
          min_dec= min(ret.corr.data.field_dec);
          max_dec= max(ret.corr.data.field_dec);
          dso=find(...
              min_ra <= self.catalogs.deep_sky_objects.RA ...
            &           self.catalogs.deep_sky_objects.RA  <= max_ra ...
            & min_dec<= self.catalogs.deep_sky_objects.DEC ...
            &           self.catalogs.deep_sky_objects.DEC <= max_dec);
          for index=1:numel(dso)
            % [ X Y ] = ad2xy*[ RA DEC ]
            st = dso(index);
            ra = self.catalogs.deep_sky_objects.RA(st);
            dec= self.catalogs.deep_sky_objects.DEC(st);
            
            xy = self.result.ad2xy*[ ra dec ]'; x = xy(1); y = xy(2);
            %[x,y] = compute_ra2xy(ra, dec, self.result.wcs.meta);
            
            [ra_h, ra_min, ra_s]      = getra(ra/15);
            [dec_deg, dec_min, dec_s] = getdec(dec);
            name = self.catalogs.deep_sky_objects.NAME{st};
            typ  = self.catalogs.deep_sky_objects.TYPE{st};
            mag  = self.catalogs.deep_sky_objects.MAG(st);
            sz   = self.catalogs.deep_sky_objects.SIZE(st); % arcmin
            
            if isfinite(sz) && sz > 1
              h = plot(x,y,'co'); 
              set(h, 'MarkerSize', ceil(sz));
            else
              h = plot(x,y,'bs'); 
            end
            hcmenu = uicontextmenu;
            
            uimenu(hcmenu, 'Label', [ 'RA=  ' num2str(ra_h)    ':' num2str(ra_min)  ':' num2str(ra_s) ]);
            uimenu(hcmenu, 'Label', [ 'DEC= ' num2str(dec_deg) ':' num2str(dec_min) ':' num2str(dec_s) ]);
            uimenu(hcmenu, 'Label', [ '<html><b>' name '</html></b>' ], 'Separator','on');
            uimenu(hcmenu, 'Label', [ 'TYPE: ' typ ]);
            if isfinite(mag)
              uimenu(hcmenu, 'Label', [ 'MAGNITUDE= ' num2str(mag)  ]);
            end
            set(h, 'UIContextMenu', hcmenu);
            
          end
        end
      end
    end % image
    
    
  end % methods
  
end % astrometry

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

function [ra_h, ra_min, ra_s] = getra(ra)
  % getra: convert any input RA (in hours) into h and min
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
  else
    disp([ mfilename ': invalid RA.' ])
    disp(ra)
  end
  if nargout == 1
    ra_h = ra_h+ra_min/60 + ra_s/3600;
  elseif nargout == 2
    ra_min = ra_min + ra_s/60;
  elseif nargout == 3
    ra_s   = abs(ra_min - fix(ra_min))*60;
    ra_min = fix(ra_min);
  end
end % getra

function [dec_deg, dec_min, dec_s] = getdec(dec)
  % getdec: convert any input DEC into deg and min
  if ischar(dec)
    dec = repradec(dec);
  end
  dec_s   = 0;
  if isstruct(dec) && isfield(dec, 'DEC')
    dec = dec.DEC;
  elseif isstruct(dec) && isfield(dec, 'deg') && isfield(dec, 'min')
    dec_deg = dec.deg;
    dec_min = dec.min;
  end
  if isnumeric(dec)
    if isscalar(dec)
      dec_deg = floor(dec); dec_min = abs(dec - dec_deg)*60;
    elseif numel(dec) == 2
      dec_deg = dec(1);   dec_min = abs(dec(2));
    elseif numel(dec) == 3
      dec_deg = dec(1);   dec_min = abs(dec(2))+abs(dec(3)/60);
    end
  else
    disp([ mfilename ': invalid DEC' ])
    disp(dec)
  end
  if nargout == 1
    dec_deg = dec_deg + dec_min/60 + dec_s/3600;
  elseif nargout == 2
    dec_min = dec_min + dec_s/60;
  elseif nargout == 3
    dec_s   = abs(dec_min - fix(dec_min))*60;
    dec_min = fix(dec_min);
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

function ret = getresult(d)
  % getresult: extract WCS and star matching information from the output files.
  %
  % input:
  %   d: directory where astrometry.net results are stored.
  
  if isdeployed || ~usejava('jvm') || ~usejava('desktop')
    disp([ '  Results are in: ' d ]);
  else
    disp([ '  Results are in: <a href="' d '">' d '</a>' ]);
  end

  if exist(fullfile(d, 'results.wcs'), 'file')
    ret.wcs  = read_fits(fullfile(d, 'results.wcs'));
    % get image center and print it
    if isfield(ret.wcs,'meta') && isfield(ret.wcs.meta,'CRVAL1')
      ret.ra  = ret.wcs.meta.CRVAL1;
      [ret.ra_h, ret.ra_min, ret.ra_s] = getra(ret.ra/15);
      ret.dec  = ret.wcs.meta.CRVAL2;
      [ret.dec_deg, ret.dec_min, ret.dec_s] = getdec(ret.dec);
      disp(  'Field center:')
      disp([ '  RA=  ' num2str(ret.ra_h)    ':' num2str(ret.ra_min)  ':' num2str(ret.ra_s) ])
      disp([ '  DEC= ' num2str(ret.dec_deg) ':' num2str(ret.dec_min) ':' num2str(ret.dec_s) ])
      % compute pixel scale
      ret.pixel_scale = sqrt(abs(ret.wcs.meta.CD1_1 * ret.wcs.meta.CD2_2 - ret.wcs.meta.CD1_2 * ret.wcs.meta.CD2_1))*3600; % in arcsec/pixel
      disp([ 'Pixel scale: ' num2str(ret.pixel_scale) ' [arcsec/pixel]' ]);
    end
  end
  if exist(fullfile(d, 'results.rdls'), 'file')
    ret.rdls = read_fits(fullfile(d, 'results.rdls'));
  end
  if exist(fullfile(d, 'results.corr'), 'file')
    ret.corr = read_fits(fullfile(d, 'results.corr'));
    % compute the transformation matrix
    if isfield(ret.corr.data, 'field_x')
      x   = ret.corr.data.field_x;
      y   = ret.corr.data.field_y;
      ra  = ret.corr.data.field_ra;
      dec = ret.corr.data.field_dec;
      AD = [ ra dec ]'; % entries (stars) are columns. 2 rows [ A ; D ]
      XY = [ x  y ]';   % entries (stars) are columns. 2 rows [ X ; Y ]
      ret.ad2xy = XY*pinv(AD);  % [ X ; Y ] = ad2xy*[ RA ; DEC ]
    end
  end
  if exist(fullfile(d, 'results.json'), 'file')
    ret.json = loadjson(fullfile(d, 'results.json'));
  end
  ret.dir     = d;
  
end % getresult

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

function [X,Y] = compute_ra2xy(RA, Dec, wcs)
% compute_ra2xy(ra, dec, wcs)
% 
% input:
%   ra,dec: coordinates to convert (in deg)
%   wcs:    WCS structure with e.g. fields CRVAL1, CRVAL2, CRPIX1, CRPIX2
%             CD1_1, CD1_2, CD2_1, CD2_2

% use: http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/wirwolf/docs/CD_PV_keywords.pdf

  CRVAL1 = wcs.CRVAL1;
  CRVAL2 = wcs.CRVAL2;
  
  CRPIX1 = wcs.CRPIX1;
  CRPIX2 = wcs.CRPIX2;
  
  CD1_1  = wcs.CD1_1;
  CD1_2  = wcs.CD1_2;
  CD2_1  = wcs.CD2_1;
  CD2_2  = wcs.CD2_2;
  
  % first convert (RA,DEC) to (eta,xi). Eq (A32-33)
  eta = (1 - tand(CRVAL2)*cosd(RA - CRVAL1)/ tand(Dec) ) ...
      / (    tand(CRVAL2)+cosd(RA - CRVAL1)/ tand(Dec) ) ;
  xi  =  tand(RA - CRVAL1)*cosd(CRVAL2)*(1 - eta*tand(CRVAL2));
  
  % then compute (X,Y)
  % use the 'no distorsion' case
  x = xi; y = eta;
  
  CD = [ CD1_1 CD1_2 ; CD2_1 CD2_2 ];
  DC = inv(CD);
  
  DC1_1 = DC(1,1);
  DC1_2 = DC(1,2);
  DC2_1 = DC(2,1);
  DC2_2 = DC(2,2);
  
  X = DC1_1*x + DC1_2*y + CRPIX1;
  Y = DC2_1*x + DC2_2*y + CRPIX2;
  
end % compute_ra2xy

function [RA,Dec] = compute_xy2ra(x,y, wcs)
% compute_xy2ra(x, y, wcs)
%
% input:
%   x,y:    pixel coordinates to convert
%   wcs:    WCS structure with e.g. fields CRVAL1, CRVAL2, CRPIX1, CRPIX2
%             CD1_1, CD1_2, CD2_1, CD2_2

% use: http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/wirwolf/docs/CD_PV_keywords.pdf

  CRVAL1 = wcs.CRVAL1;
  CRVAL2 = wcs.CRVAL2;
  
  CRPIX1 = wcs.CRPIX1;
  CRPIX2 = wcs.CRPIX2;
  
  CD1_1  = wcs.CD1_1;
  CD1_2  = wcs.CD1_2;
  CD2_1  = wcs.CD2_1;
  CD2_2  = wcs.CD2_2;
  
  x = CD1_1*(X - CRPIX1) + CD1_2*(Y - CRPIX2);
  y = CD1_1*(X - CRPIX1) + CD2_2*(Y - CRPIX2);
  
  % use the 'no distorsion' case
  xi = x; eta = y;
  
  % Eq (A29-31)
  alpha = atand( xi/cosd(CRVAL2)/(1-eta*tand(CRVAL2)) );
  RA    = alpha+CRVAL1;
  Dec   = atand( (eta+tand(CRVAL2))*cosd(alpha) / ( 1 - eta*tand(CRVAL2) )  );
  
end % compute_xy2ra
