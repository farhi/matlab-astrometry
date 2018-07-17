function [s,f,m] = peakwidth(im, x0, dx)
% peakwidth: determine a peak width along X and Y
%
% input:
%   im:     image (rgb)
%   x0:     peak location guess, or empty to use max
%   dx:     pixel number to extract around peak (optional, default=+/-10)
% output:
%   s:      width    along dimensions (vector)
%   f:      centroid along dimensions (vector)

  % when x,y are given we search around this location within dx
  % else we search the max in image
  im1=im; f = 0; s = 0; m=0;
  if ndims(im1) == 3, im1 = rgb2gray(im1); end
  if nargin < 3,  dx=''; end
  if isempty(dx), dx=20; end
  if nargin < 2,  x0=''; end
  if numel(x0) < 2
    % when location is not given we use the max of whole image
    [x0,m] = peak_max(im1);
  else
    % extract sub-image around x0 within dx, then get the max location
    X=(x0(1)-dx):(x0(1)+dx); X=round(X(X>=1 & X<=size(im1,1)));
    Y=(x0(2)-dx):(x0(2)+dx); Y=round(Y(Y>=1 & Y<=size(im1,2)));
    if isempty(X) || isempty(Y), return; end
    [x0,m] = peak_max(im1(X,Y),'proj');
    x0 = x0 + [ min(X) min(Y) ] -1;
  end
  % extract the portion of the image to analyze around x0
  if numel(x0) < 2, return; end
  dx = ceil(dx);
  X  = ceil((x0(1)-dx):(x0(1)+dx)); 
  Y  = ceil((x0(2)-dx):(x0(2)+dx)); 
  % restrict to usable range
  X  = X(X>=1 & X<=size(im1,1));
  Y  = Y(Y>=1 & Y<=size(im1,2));
  if isempty(X) || isempty(Y), return; end
  
  % determine centroid and width (gaussian)
  im2=double(im1(X,Y));
  [s1, f1] = width1(X, im1(X,     x0(2)));
  [s2, f2] = width1(Y, im1(x0(1), Y));
  
  is = find(im2-min(im2(:)) > max(im2(:))/2); % nb of pixels above half height
  s  = sqrt(numel(is)); 
  s  = [ s s ];
  f  = [ f1 f2 ];
  if im1(x0(1), x0(2))*.8 > im1(round(f(1)), round(f(2))), f = x0; end
  % s = [ s1 s2 ];
  
end % peakwidth

% ------------------------------------------------------------------------------

function [centroid, int] = peak_max(im,method)
  % peak_max: get max location in image
  %
  % input:
  %  im:     image (MxN matrix luminance)
  %  method: 'proj' or 'max' (default)
  %
  % output:
  %  centroid: index location of max (vector)
  %  int:      intensity at max
  
  if nargin < 2, method = ''; end
  
  if strcmp(method, 'proj')
    % we use the summed image (less noisy)
    [int, x1] = max(sum(im, 2));
    [int, x2] = max(sum(im, 1));
  else
    [int,x1]  = max(im(:));
    [x1,x2]   = ind2sub(size(im),x1);
  end
  centroid  = [ x1 x2 ];
  int       = double(int);
end % 
