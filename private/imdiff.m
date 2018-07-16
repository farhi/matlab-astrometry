function [ret_t, ret_R, theta] = imdiff(self, img1, img2)
      % diff: compute the difference between img2 and img1 used as reference.
      %
      %   [ret_t, theta] = diff(self, img1, img2)
      %     returns the translation and rotation between img1 and img2.
      %   [ret_t, theta] = diff(self, img2)
      %     use the reference 'light' image as img1.
      
      ret_t=[]; ret_R = []; theta=[];
      if nargin < 3, return; end
      if isempty(img1) || isempty(img2), return; end
      
      % the images should have the same scaling
      tol_trans = self.toleranceTranslation;      % in percent
      tol_trans = tol_trans*min(img2.image_size); % in pixels
      tol_trans = max(tol_trans, 20);
      tol_rot   = self.toleranceRotation;         % in deg
      [x1,y1,x2,y2, p1_orig,p2_orig, p1_axis,p2_axis] = ...
        analyse_dist_angles(img1.points, img2.points, tol_trans, tol_rot);
      
      % get the best similarity
      if ~p1_orig || ~p2_orig || ~p1_axis || ~p2_axis
        disp([ mfilename ': WARNING: not enough common control points for ' img1.id ' ' img2.id ' (axis).' ])
        return
      end
  
      % identify the current points which are common to the reference, and match
      % distances AND rotations
      p1 = p1_orig; p2=p2_orig;
      d1  = sqrt( (x1-x1(p1)).^2 + (y1-y1(p1)).^2 );
      t1  = atan2( y1-y1(p1),       x1-x1(p1))*180/pi;  % deg
      d2  = sqrt( (x2-x2(p2)).^2 + (y2-y2(p2)).^2 );
      t2  = atan2( y2-y2(p2),       x2-x2(p2))*180/pi;  % deg
      [ok1,ok2] = find_similar2(d1,             d2,             tol_trans, ...
                                t1-t1(p1_axis), t2-t2(p2_axis), tol_rot);
                                
      if numel(ok1) <= 1 || numel(ok2) <= 1
        disp([ mfilename ': WARNING: dist/angle mismatch between control points for ' img1.id ' ' img2.id ' (common).' ])
        return
      end
  
      % we make a check for wrong guesses
      delta = (t1(ok1)-t1(p1_axis)) - (t2(ok2)-t2(p2_axis));
      bad   = find(abs(delta) > tol_rot);
      ok1(bad) = [];
      ok2(bad) = [];
      if numel(ok1) <= 1 || numel(ok2) <= 1
        disp([ mfilename ': WARNING: angle mismatch between control points for ' img1.id ' ' img2.id ' (common2). Had initially ' num2str(numel(delta)) ' points.'])
        return
      end
  
      % compute the mean translation (x,y) and rotation angle wrt to reference image
      x1  = x1(ok1); y1=y1(ok1);
      x2  = x2(ok2); y2=y2(ok2);
      
      % theta2-theta1 is an estimate of the rotation angle
      % compute the affine transformatin
      [ret_R, ret_t] = rigid_transform_3D([x2 ; y2]', [x1 ; y1]');
      if isempty(ret_R)
        disp([ mfilename ': WARNING: invalid affine transformation. Skipping.'])
        ret_t = [];
        return
      end
      % compute an estimate of the translation and rotation from the identified
      % control points orig and axis
      theta1 = t1(p1_axis); % t1(p1_orig) == 0
      theta2 = t2(p2_axis);
      theta = atan2(ret_R(2,1),ret_R(1,1))*180/pi;
      if theta < -90 && theta2-theta1 > 90
        theta = theta + 360;
      elseif theta > 90 && theta2-theta1 < -90
        theta = theta - 360;
      end
      if abs(theta-(theta2-theta1)) > tol_rot
        disp([ mfilename ': WARNING: invalid affine rotation theta=' num2str([theta theta2-theta1]) '. Skipping.']);
        ret_t = [];
        return
      end

    end % diff
