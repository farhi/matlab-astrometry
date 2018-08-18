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
  
    function astrometry(filename, varargin)
      % astrometry: loads an image and identifies its objects using astrometry.net
      % 
      self.executables = find_executables;
      
    end % astrometry
    
    function urlread(self, filename, varargin)
      %  input:
      %    filename: an image to annotate
      
      % varargin: ra, dec, radius, scale-lower, scale-upper
      
      if ismac,      precmd = 'DYLD_LIBRARY_PATH= ;';
      elseif isunix, precmd = 'LD_LIBRARY_PATH= ; '; 
      else           precmd=''; end
      
      % python client.py -k API_KEY -u IMAGE_FILE -w --newfits=results.fits --kmz=results.kmz --wcs=results.wcs -a results.json

    
    end % urlread
    
    function system(self, filename, varargin)
      %  input:
      %    filename: an image to annotate
      %    any name/value pairs as:
      %      ra:      approximate RA coordinate  (e.g. deg or 'hh:mm:ss')
      %      dec:     approximate DEC coordinate (e.g. deg or 'deg:mm:ss')
      %      radius:  approximate field size      (in deg)
      %      kmz:     output KML file             (requires wcs2kml)
      %      depth:   number of field objects to look at (e.g. 1-50)
      
      % varargin: ra, dec, radius, scale-low, scale-high
      
      if ismac,      precmd = 'DYLD_LIBRARY_PATH= ;';
      elseif isunix, precmd = 'LD_LIBRARY_PATH= ; '; 
      else           precmd=''; end
      
      % --dir tempdir
      % solve-field IC5146-Cocon-2018-08-11.jpg  --ra 21:53:24 --dec 47:16:14 --radius 2  -l 600 -L 0.5 -H 2 --overwrite
    
    end % system
    
    
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
  for exe =  { 'solve-field','sextractor','python', 'wcs2kml', 'client.py','googleearth','google-earth-pro' }
  
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

  this_path = fullfile(fileparts(which(mfilename)));
