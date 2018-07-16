function [Y,M] = imaffine(X, ret_R, ret_t)
  % apply an affine transform to an image
  %
  % input:
  %  im: initial image
  %  R,t: rotation matrix(2x2) and translation
  %
  % output:
  %  Y: transformed image
  %  M: mask indicating assigned pixels in new image
	m = size(X,1);
	n = size(X,2);
	s = size(X); if numel(s) < 3, s(3)=1; end
	x1 = 1:m;
	y1 = 1:n;
	[x1,y1]= meshgrid(x1,y1);
	A = [ x1(:) y1(:) ];
	B = (ret_R*A') + repmat(ret_t, 1, size(A,1));
	clear A
	B = round(B); x2=B(1,:); y2=B(2,:); 
	clear B
	Y = X*0;
	M=uint8(zeros(m,n));
	x2(x2<1)= 1; y2(y2<1)=1;
	x2(x2>size(X,1)) =1;
	y2(y2>size(X,2)) =1;
	for z=1:s(3)
	  i1 = sub2ind(size(X),x1,y1,z*ones(size(x1)));
	  i2 = sub2ind(size(Y),x2,y2,z*ones(size(x2)));
	  Y(i2) = X(i1);
	  if z == 1, M(i2) = 1; end
	end
end % imaffine
