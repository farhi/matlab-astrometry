function points = find_control_points(im, N, deadPixelArea)
  % find_control_points: find N points separated with tol_trans*10 pixels
  if nargin < 2, N=30; end
  if ndims(im) == 3, im = rgb2gray(im); end
  
  points = [];
  if isempty(im), return; end
  
  points.x  = [];
  points.y  = [];
  points.m  = [];
  points.sx = [];
  points.sy = [];
  
  % we need N points within the image. Area per point is prod(size(im))/N
  % exclusion distance is sqrt(prod(size(im))/N) (full width)
  dx = sqrt(prod(size(im))/N)/2;  % and we half it

  % find 'max' intensity locations. Every time we have a spot, we 'blank' it
  % around so that other ones are separated by 'tol_trans*10'
  for p=1:N
    [s,f,m, im, iter] = max_and_zero(im, dx/2, deadPixelArea);
    if isempty(f) || any(f == 0), continue; end % can not find a new star
    if numel(f) == 2 && numel(points.x) && points.x(end) == f(1) && points.y(end) == f(2)
      continue;
    end
    points.x(end+1) = f(1);
    points.y(end+1) = f(2);
    points.m(end+1) = m;
    points.sx(end+1)= s(1);
    points.sy(end+1)= s(2);
  end

  % sort points with increasing 'y' value
  [points.y,p] = sort(points.y);
  points.x = points.x(p);
  points.m = points.m(p);
  points.sx= points.sx(p);
  points.sy= points.sy(p);
  points.handle = [];
end % find_control_points
