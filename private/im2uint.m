function im = im2uint(im, cl)
  % im2uint: cast the image to uint class
  
  if nargin < 2, cl='uint8'; end
  if isempty(im), return; end
  if strcmp(class(im), cl), return; end
  
  switch cl
  case 'uint8'
    cl = 8;
  case 'uint16'
    cl = 16;
  case 'single'
    im = single(imdouble(im));
    return
  case 'double'
    im = imdouble(im);
    return
  otherwise
    error([ mfilename ': can only cast image into uint8 and uint16, not ' num2str(cl) ]);
  end
  
  im = cast(imdouble(im)*2^cl, [ 'uint' num2str(cl) ]);
      
