function [s,f] = width1(x,s)
% width1: compute the width and central axis value of a 1D signal
%
% input:
%   x: axis (vector)
%   s: signal (vector)
%
% output:
%   s: second moment (width)
%   f: first moment (mean axis value)
  x=x(:); s=double(s(:)); s=s-min(s(:));
  s=s/sum(s);
  
  % first moment (mean)
  f = sum(s.*x); % mean value
  
  % second moment: sqrt(sum(x^2*s)/sum(s)-fmon_x*fmon_x);
  s = sqrt(sum(x.*x.*s) - f*f);
end % width1
