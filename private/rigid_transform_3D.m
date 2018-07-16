function [R,t] = rigid_transform_3D(A, B)
  % http://nghiaho.com/?page_id=671
  % apply affine transform to coordinates with:
  %   B = (ret_R*A') + repmat(ret_t, 1, size(A,1))
    if nargin ~= 2
	    error('Missing parameters');
    end

    if ~(all(size(A) == size(B)))
      R=[]; t=[];
      return;
    end

    centroid_A = mean(A);
    centroid_B = mean(B);

    N = size(A,1);

    H = (A - repmat(centroid_A, N, 1))' * (B - repmat(centroid_B, N, 1));

    [U,S,V] = svd(H);

    R = V*U';

    if det(R) < 0 && size(V,2) > 3
        disp('Reflection detected');
        V(:,3) = -V(:,3); 
        R = V*U';
    end

    t = -R*centroid_A' + centroid_B';
  end % rigid_transform_3D
