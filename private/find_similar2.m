function [index_A_in_B, index_B_in_A] = find_similar2(A, B, find_tol1, C, D, find_tol2)
  % find_similar2: find indices in 'A' and 'C' that are similar to 'B' and 'D'
  %   within 'find_tol'
  %
  %   A(index_A_in_B) ~ B
  %   B(index_B_in_A) ~ A
  %
  % input:
  %   A:        vector for which values must be compared with the reference
  %   B:        list of 'reference' values
  %   find_tol1: tolerance used for comparison between A and B
  %   C:        vector for which values must be compared with the reference
  %   D:        list of 'reference' values
  %   find_tol2: tolerance used for comparison between C and D
  %
  % output:
  %   index_A_in_B: 'A' elements which are in 'B' within 'find_tol'
  %   index_B_in_A: 'B' elements which are in 'A' within 'find_tol'

  iA = zeros(size(A)); iB = zeros(size(B));
  index_A_in_B = [];
  index_B_in_A = [];
  
  for indexB=1:numel(B)
    x = B(indexB);
    [mindiff1, fA] = min(abs(A-x));
    y = D(indexB);
    [mindiff2, fC] = min(abs(C-y));
    if mindiff1 < find_tol1 && mindiff2 < find_tol2 && ~iA(fA) && ~iB(indexB)
      iA(fA) = 1; iB(indexB)=1;
      index_A_in_B = [ index_A_in_B fA ];
      index_B_in_A = [ index_B_in_A indexB ];
    end
  end
end % find_similar2

