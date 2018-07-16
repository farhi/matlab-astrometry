function [im, img, exif] = imread_single(source, images, flag)
  % imread_single: read a single image.
  %
  %   source: the image filename, or matrix, or image index
  %   images: the loaded images (structure array)
  %   flag:   when 0, the image is not read when not new
  
  if nargin < 3, flag=1; end

  im = []; img = []; exif = [];
  if iscellstr(source), source = char(source); end
  if isempty(source), return; end
  % new image
  if ischar(source) && strncmp(source,'file://',7)
    source = source(8:end);
  end
  if ischar(source) && exist(source, 'file')
    % check if image exists
    try
      exif     = imfinfo(source);
      if flag, im       = imread(source); end
      disp([ mfilename ': Reading ' source ' (from file)']);
    catch ME
      disp([ mfilename ': Skipping ' source ]);
      return
    end
  elseif (isnumeric(source) || isstruct(source)) && isscalar(source)
    % source is an image index or struct
    if isnumeric(source) && source && numel(images) >= source
      img    = images(source);
      source = img.source;
    elseif isstruct(source)
      img    = source;
      source = img.source;
    else img = [] ; end
    % make sure we can get the image (matrix), when img was retrieved but not image
    if ~isempty(img) && flag && ischar(img.source)
      disp([ mfilename ': Reading ' img.source ' (from image ref)' ]);
      im = imread(img.source);
    end
  elseif isnumeric(source) && ~isempty(source)
    exif     = [];
    im       = source;
  else
    if ischar(source) || isscalar(source)
      source
    end
    error([ mfilename ': argument should be char or matrix. Is ' class(source) ' [' num2str(numel(source)) ']' ])
    
    return
  end

  if flag && isempty(im) return; end

  % create fake exif if not set
  if isempty(exif)
    if ischar(source)
      exif.Filename     = source;
    else
      exif.Filename     = sprintf('Image_%d', numel(images)+1);
    end
    exif.Format       = 'jpg';
    exif.FormatVersion= '';
    exif.FileModDate  = datestr(now);
    exif.Software     = mfilename;
    
    exif.FileSize     = numel(im);
    exif.Width        = size(im,2);
    exif.Height       = size(im,1);
    exif.BitDepth     = size(im,3)*8;
    if size(im,3) == 1
      exif.ColorType    = 'grayscale';
    else
      exif.ColorType    = 'truecolor';
    end
  end
  if numel(exif) > 1, exif = exif(1); end
  if isfield(exif, 'DigitalCamera') && isfield(exif.DigitalCamera, 'ExposureTime')
    exif.ExposureTime = exif.DigitalCamera.ExposureTime;
  elseif isfield(exif, 'Shutter')
    exif.ExposureTime = str2num(strtok(exif.Shutter));
  end
  
  if ~isstruct(img)
    % create the img structure holding information and EXIF
    img.index       = numel(images)+1;
    [~,img.id]      = fileparts(exif.Filename);
    img.source      = source;
    img.image       = []; % we do not store the images. Many will be huge.
    img.exif        = exif;
    img.width       = 0;
    img.image_size  = 0;
    img.image_sum   = 0;
    img.sharpness   = 0;
    img.intensity   = 0;
    img.type        = [];
    img.points      = struct('x',[],'y',[],'handle',[],'m',[], ...
      'sx',[],'sy',[], 'sharpness', []);
  end
  if ~isempty(im)
    img.image_size  = size(im);
    img.image_sum   = sum(double(im(:)));
  end
  if all(img.image_size == 0) && isfield(img.exif,'Width') && isfield(img.exif,'Height')
    img.image_size = [ img.exif.Width img.exif.Height ];
  end
  
end % imread_single
