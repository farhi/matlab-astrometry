classdef astrometry < handle
  % astrometry: A Matlab class to annotate astrophotography images (identify objects/astrometry)
  %
  %  Purpose
  %
  % (c) E. Farhi, 2018. GPL2.

  properties
    toleranceTranslation  = 0.01; % in percent
    toleranceRotation     = 1;    % in deg
    deadPixelArea         = 9;
    
    images = [];
    %      id
    %      source (or matrix if directly given)
    %      image_size
    %      image_sum
    %      exif
    %      type (light, dark, flat, skip, none)
    %      points
    %      rotation
    %      translation
    %      sharpness
    %      width
    %      thumbnail
    
    nbControlPoints = 50;
    UserData    = [];
    sensor_size = [ 25.1 16.7 ];
    focal_length= 1200;
    
    catalogs = [];
  end
  
  properties (Access=private)
    figure       = [];
    
  end % properties
  
  methods
  
    % I/O stuff ================================================================
    function self = astrometry(source)
      % astrometry: create the astrometry session, load file
      %
      % self = astrometry(image)
      %   import given image, specified as filenames or matrix (supports wildcards)

      % load catalogs: objects, stars for astrometry
      disp([ mfilename ': Welcome !' ]);
      
      if nargin == 1 % import images
        img = load(self, source);
      end
      
      self.catalogs = load(mfilename);
      
      % display available catalogs
      for f=fieldnames(self.catalogs)'
        name = f{1};
        if ~isempty(self.catalogs.(name))
          num  = numel(self.catalogs.(name).RA);
          if isfield(self.catalogs.(name), 'Description')
            desc = self.catalogs.(name).Description;
          else desc = ''; end
          disp([ mfilename ': ' name ' with ' num2str(num) ' entries.' ]);
          disp([ '  ' desc ])
        end
      end
      
    end % astrometry
    
    function [im, img] = imread(self, source, flag)
      % imread: read an image, and store its information
      %
      % Supported file formats include JPG, PNG, BMP, GIF, TIFF, FITS
      %
      %  [im, img] = imread(self, filename or matrix or index)
      %    returns the image matrix and information
      %  [~, img] = imread(self, filename or matrix or index, 0)
      %    only returns the image information, avoiding to read the file
      %
      % see also: imread
      
      % flag: when 0, do not force to read matrix images
      if nargin < 3, flag=1; end
      
      % first check if image(s) are already loaded
      [found, src] = exist(self, source);
      if ~iscell(src) && ~isstruct(src), src = { src }; end
      img   = []; im = {};
      if isempty(src), return; end
      
      for index=1:numel(found)
      
        % read new image, or retrieve it
        if isstruct(src)
          this_src = src(index);
        else
          this_src = src{index};
        end
        [this_im, this_img] = imread_single(this_src, self.images, flag);
        if ischar(this_src)
          if strfind(lower(this_src), 'dark'), this_img.type = 'dark';
          elseif strfind(lower(this_src), 'flat'), this_img.type = 'flat';
          end
        end
        
        if isempty(this_img) continue; end  % invalid image
        
        % set image index
        if isnan(found(index))
          this_img.index = numel(self.images)+1;
        else % update
          this_img.index = found(index);
        end
        % store
        if isempty(self.images) && this_img.index > 0
          self.images    = this_img;
        elseif this_img.index > 0
          self.images(this_img.index) = this_img; 
        else
          return
        end

        % add to return arguments
        if isempty(img), img        = this_img;
        else             img(end+1) = this_img; end
        if flag, im{end+1} = this_im; end
      end % for
      
      if iscell(im) && numel(im) == 1
        im = im{1};
      end
      
    end % imread
    
    function img = load(self, source)
      % load: read an image, and store its information
      %
      % Supported file formats include JPG, PNG, BMP, GIF, TIFF, FITS
      %
      %  img = load(self, filename or matrix or index)
      %    returns the image matrix and information.
      [~,img] = imread(self, source, 0);
    end % load
    
    function [found, src] = exist(self, source)
      % exist: check if an image is already loaded
      %
      %   found = exist(self, source)
      %     returns the index of any matching image, or nan (not found)

      % handle array of images
      found = []; src = {};
      if ischar(source)
        source = resolvefiles(source);
      end
      if iscellstr(source) && numel(source) == 1
        source = source{1}; 
      end
      if (iscell(source) && numel(source) > 1) || (isnumeric(source) && ~isscalar(source) ...
        && numel(source) == max(size(source))) || (isstruct(source) && numel(source) > 1)
        for index=1:numel(source)
          if iscell(source), this=source{index};
          else this=source(index); end
          [found(end+1), src{end+1}] = exist(self, this);
        end
        return
      end

      found = nan; src = source;
      if isnumeric(source) && isscalar(source) && source <= numel(self.images)
        found = source;
      else
        for index=1:numel(self.images)
          this = self.images(index);  % a struct
          if isempty(this), continue; end
          if ischar(source) && ~isempty(this.source) && ischar(this.source)
            % found: same name
            if strcmp(this.source, source) || strcmp(this.id, source)
              found = index;
              break; 
            end
          elseif isnumeric(source)
            % found: equal sum and size
            if all(size(this.image_size) == size(source)) ...
              && this.image_sum && sum(source(:)) == this.image_sum
              found = index;
              break; 
            end
          elseif isstruct(source)
            if any(strcmp(this.source, { source.source })) ...
            || any(strcmp(this.id, { source.id }))
              found = index;
              break; 
            end
          end
        end % for
      end
      
    end % exist
   
    
    % compute stuff ============================================================
   
    
    function [ret_t, ret_R, self] = diff(self, img1)
      % diff: compute the translation and rotation wrt reference image
      %
      % [t,R] = diff(self, image)
      %   return the translation and rotation of image wrt reference
      [~,img1] = imread(self, img1, 0);
      [~,img2] = imread(self, self.reference, 0);
      ret_t= []; ret_R = [];
      if isempty(img1) || isempty(img2), return; end
      if numel(img1.points.x) < 2
        [img1, self] = cpselect(self, img1);
      end
      if numel(img2.points.x) < 2
        [img2, self] = cpselect(self, img2);
      end
      
      if isempty(img2), return; end
      [ret_t, ret_R] = imdiff(self, img2, img1);
    end % diff
    
  
    
    function [this_img, self] = align(self, img, im)
      % align: automatically set control points
      [this_img, self] = cpselect(self, img, im);
    end
    
    function [nimg, self] = cpselect(self, img, im)
      % cpselect: automatically set control points
      %
      % cpselect(self)
      %   set all control points and compute sharpness
      % cpselect(self, images)
      %   set control points and compute sharpness on given images
      %   'images' can be given as name or index
      
      if nargin < 2, img=1:numel(self.images); end
      
      if numel(img) > 1
        delete(findall(0, 'Tag', [ mfilename '_waitbar' ]));
        wb  = waitbar(0, [ mfilename ': Aligning images (close to abort)...' ]); 
        set(wb, 'Tag', [ mfilename '_waitbar' ]);
        disp([ mfilename ':   Aligning... ' datestr(now) ])
      else wb = [];
      end
      t0 = clock; 
      
      nimg = [];
      
      for index=1:numel(img)
        this_img = img(index);
        
        % get the image
        [~, this_img]  = imread(self, this_img, 0);
        if isempty(this_img.points.x) || nargin < 2
        
          if nargin < 3 || isempty(im) || ~isnumeric(im)
            [im, this_img] = imread(self, this_img, 1);
          end
          im = rgb2gray(im);
          
          % update waitbar and ETA display
          if numel(img) > 1 && ~isempty(wb)
            if ~ishandle(wb)
              disp('Align: Aborting (user closed the waitbar).')
              break;
            end
            % compute ETA
            dt_from0     = etime(clock, t0);
            dt_per_image = dt_from0/index;
            % remaining images: numel(s)-index
            eta    = dt_per_image*(numel(img)-index+1);
            ending = addtodate(now, ceil(eta), 'second');
            ending = [ 'Ending ' datestr(ending) ];
            eta    = sprintf('ETA %i [s]. %s', round(eta), ending);
            waitbar(index/numel(img), wb, [ 'Align: ' this_img.id ' (close to abort)...' ]);
            try
            set(wb, 'Name', [ num2str(index) '/' num2str(numel(img)) ' ' eta ]);
            end
          end
          
          this_img.points = ...
              find_control_points(im, self.nbControlPoints, self.deadPixelArea);
          this_img.width  = ...
              sqrt(sum(this_img.points.sx.^2.*this_img.points.sy.^2)) ...
              /numel(this_img.points.sx);
          this_img.sharpness  = ...
              sqrt(sum(this_img.points.sharpness.^2)) ...
              /numel(this_img.points.sharpness);
          this_img.intensity  = ...
              sum(this_img.points.m) ...
              /numel(this_img.points.m);
        end
        
        if isempty(nimg), nimg = this_img;
        else nimg = [ nimg this_img]; end
        self.images(this_img.index) = this_img;
      end % for
      if numel(img) > 1 && ~isempty(wb)
        delete(findall(0, 'Tag', [ mfilename '_waitbar' ]));
        disp([ mfilename ':   Alignment done ' datestr(now) ])
      end
      
    end % cpselect
    
    function solve(self, img, ra, dec, radius)
      % astrometry: determine the location of an image in RA/DEC coordinates
      
      % need focal and sensor dimensions OR field-of-view
      % need an estimate of the RA/DEC centre, and search radius
      if nargin < 3, ra=[]; end
      if nargin < 4, dec=[]; end
      if nargin < 5, radius=[]; end
      if ischar(ra)
        found = findobj(self, ra);
        if isempty(found)
          disp([ mfilename ': Did not find object ' ra ])
          return; 
        end
        ra = found.RA;
        dec= found.DEC;
      end
      
      if isempty(ra), ra=180; end
      if isempty(dec), dec=0; end
      if isempty(radius), radius=180; end
      
      [im,img] = imread(self, img);  % make sure we have that image. Can be a filename

      % make sure we have control points (cpselect)
      if isempty(img.points) || ~isstruct(img.points) || ...
        ~isfield(img.points, 'x') || numel(img.points.x) < 2
        % compute control points
        disp([ mfilename ':   computing control points for ' img.id ])
        [img, self] = cpselect(self, img, im);
      end
      
      % convert control point coordinates to degrees (origin=image corner, relative)
      
      % angular field of view [deg]
      FOV = 2*atan2(self.sensor_size/2,self.focal_length)*180/pi;
      
      x = img.points.x*FOV(1);  % now in DEG, starting from image corner=0
      y = img.points.y*FOV(2);
      
      % get a list of objects from catalogs in the search radius
      catalog = [];
      for f=fieldnames(self.catalogs)'
        
        index = find(abs(self.catalogs.(f{1}).RA - ra) < radius/2 ...
                   & abs(self.catalogs.(f{1}).DEC - dec) < radius/2);
        if ~isempty(index)
          for F={'RA','DEC','DIST','MAG','SIZE','TYPE','NAME'}
            this = self.catalogs.(f{1}).(F{1});
            if ~isfield(catalog, F{1})
              catalog.(F{1}) = this(index);
            else
              catalog.(F{1}) = [ catalog.(F{1}) ; this(index) ];
            end
          end
        end
        disp([ mfilename ': found ' num2str(numel(index)) ' objects from ' f{1} ])
      end

      % compute difference (diff)
      img1          = img;
      img1.points.x = x;
      img1.points.y = y;

      img2.points.x = catalog.RA';
      img2.points.y = catalog.DEC';
      img2.id       = 'Catalogs';

      [ret_t, ret_R] = imdiff(self, img2, img1);
    
    end % solve
    
    function found = findobj(self, name)
      % findobj(ma, name): find a given object in catalogs.
      found = [];
      for f=fieldnames(self.catalogs)'
        catalog = self.catalogs.(f{1});
        if ~isfield(catalog, 'MAG'), continue; end
        % search for name
        index = find(~cellfun(@isempty, strfind(catalog.NAME, [ ';' name ';' ])));
        if isempty(index)
        index = find(~cellfun(@isempty, strfind(catalog.NAME, [ name ';' ])));
        end
        if isempty(index)
        index = find(~cellfun(@isempty, strfind(catalog.NAME, [ ';' name ])));
        end
        if isempty(index)
        index = find(~cellfun(@isempty, strfind(catalog.NAME, [ name ])));
        end
        if ~isempty(index)
          found.index   = index(1);
          found.catalog = f{1};
          found.RA      = catalog.RA(found.index);
          found.DEC     = catalog.DEC(found.index);
          found.MAG     = catalog.MAG(found.index);
          found.TYPE    = catalog.TYPE{found.index};
          found.NAME    = catalog.NAME{found.index};
          break;
        end
      end
      if ~isempty(found)
        disp([ mfilename ': Found object ' name ' as: ' found.NAME ])
      end
    end % findobj
    
   
  
  end % methods
  
end % astrometry

