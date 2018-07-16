function [x1,y1,x2,y2, p1_orig,p2_orig, p1_axis,p2_axis] = analyse_dist_angles(points1, points2, tol_trans, tol_rot)

  if ~isstruct(points1) || ~isstruct(points2) || ...
    ~isfield(points1, 'x') || ~isfield(points2, 'x') || ...
    isempty(points1.x) || isempty(points2.x)
    x1=[]; y1=[]; x2=[]; y2=[]; 
    p1_orig=0; p2_orig=0; 
    p1_axis=0; p2_axis=0;
    return
  end
  % there is a global translation, and a minor rotation which brings some of
  % the current image back onto some of the previous one.
  x2 = points2.x; x1=points1.x;
  y2 = points2.y; y1=points1.y;
  i1 = find(x1>0 & y1>0);
  i2 = find(x2>0 & y2>0);
  x1 = x1(i1);  x2 = x2(i2); 
  y1 = y1(i1);  y2 = y2(i2);
  
  % we search for similar points, taking as reference one of the control points
  % the distances are p1_orig || ~p2_orig || ~p1_axis || ~p2_axis preserved in all images, and independent of rotations.
  ok1_best = 1; ok2_best = 1; p1_orig = 0; p2_orig = 0; p1_axis  = 0;  p2_axis = 0;
  for p1=1:numel(x1)
    % we select an initial control point in the reference image, that is used
    % as origin for the distances.
    % the reference can change when the second image moves and the overlap can
    % exclude some stars from the intersection.
    d1  = sqrt( (x1-x1(p1)).^2 + (y1-y1(p1)).^2 );
    
    for p2=1:numel(x2)
      % then we do the same on the current image
      d2  = sqrt( (x2-x2(p2)).^2 + (y2-y2(p2)).^2 );
      % d1(ok1) = d2(ok2)
      [ok1,ok2] = find_similar(d1, d2, tol_trans);
      
      if numel(ok1) > ok1_best && numel(ok2) > ok2_best 
        % store this solution as best bet, with at least 2 similar points
        ok1_best = numel(ok1); ok2_best = numel(ok2);
        p1_orig  = p1;         p2_orig  = p2;
      end
    end
  end
  
  if ~p1_orig || ~p2_orig, return; end
  
  % now we search for a second control point from which angles are measured
  
  ok1_best = 1; ok2_best = 1; 
  % compute absolute angle values using the same 'origin' as that for distances
  t1  = atan2( y1-y1(p1_orig), x1-x1(p1_orig))*180/pi;  % deg
  t2  = atan2( y2-y2(p2_orig), x2-x2(p2_orig))*180/pi;  % deg
  for p1=1:numel(x1)
    if p1 == p1_orig, continue; end
    for p2=1:numel(x2)
      if p2 == p2_orig, continue; end
      % we try combinations so that angle(axis-orig) is used as reference
      % and other stars angles are measured from that direction.
      % this is then independent of the image rotations.
      [ok1,ok2] = find_similar(t1-t1(p1), t2-t2(p2), tol_rot);
      % the match for distances and angle should coincide
      if numel(ok1) > ok1_best && numel(ok2) > ok2_best 
        ok1_best = numel(ok1); ok2_best = numel(ok2);
        p1_axis  = p1;         p2_axis  = p2;
      end
    end
  end
end % analyse_dist_angles
