function gray = rgb2gray (rgb)
% rgb2gray: transform an rgb image to grayscale

  gray = rgb;
  if isempty(rgb) || ndims(rgb) == 2, return; end
  
  gray = (0.2989 * rgb(:,:,1)) + (0.5870 * rgb(:,:,2)) + (0.1140 * rgb(:,:,3));
  gray = cast(gray, class(rgb));
    
end % rgb2gray
