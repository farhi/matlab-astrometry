function [s,f,m, im, iter, sharpness] = max_and_zero(im, dx, deadPixelArea, iter)
  % max_and_zero: search for a maximum intensity point and zero image around it
  %
  % input:
  %   im: image (gray m*n)
  %   dx:  area to zero after search
  %   deadPixelArea: area around dead pixels
  % output:
  %   x1,y1: location in image
  %   m1:    intensity
  
  if nargin <3, deadPixelSize=9; end
  if nargin <4, iter=1; end
  if ndims(im) == 3, im = rgb2gray(im); end
  sharpness = 0;
  
  if iter>5, s=[]; f=[]; m=[]; return; end
  
  % search for max intensity and width
  [s, f, m]  = peakwidth(im, [ ], dx);
  
  % get image around peak, and compute sharpness
  dx1 = round((f(1)-dx):(f(1)+dx));
  dy1 = round((f(2)-dx):(f(2)+dx));
  dx1=dx1(dx1>=1 & dx1 <=size(im,1));
  dy1=dy1(dy1>=1 & dy1 <=size(im,2));
  % blank image around max
  if prod(s) > 4
    sharpness   = image_sharpness(im(dx1,dy1));
  end
  im(dx1,dy1) = 0;
  
  % recursive call if that guess is not acceptable
  % remove dead pixels (too sharp peaks) and image edges.
  if prod(s) <= 9
    [s,f,m, im, iter, sharpness] = max_and_zero(im, dx, deadPixelArea, iter+1);
    return
  end

end % max_and_zero
