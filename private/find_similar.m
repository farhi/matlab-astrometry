function [index_A_in_B, index_B_in_A] = find_similar(A, B, find_tol)
  % find_similar: find indices in 'A' that are similar to 'B' 
  %   within 'find_tol'
  %
  %   A(index_A_in_B) ~ B
  %   B(index_B_in_A) ~ A
  %
  % input:
  %   A:        vector for which values must be compared with the reference
  %   B:        list of 'reference' values
  %   find_tol: tolerance used for comparison
  %
  % output:
  %   index_A_in_B: 'A' elements which are in 'B' within 'find_tol'
  %   index_B_in_A: 'B' elements which are in 'A' within 'find_tol'
  
  iA = zeros(size(A)); iB = zeros(size(B));
  index_A_in_B = [];
  index_B_in_A = [];
  
  for indexB=1:numel(B)
    x = B(indexB);
    [mindiff, fA] = min(abs(A-x));
    if mindiff < find_tol && ~iA(fA) && ~iB(indexB)
      iA(fA) = 1; iB(indexB)=1;
      index_A_in_B = [ index_A_in_B fA ];
      index_B_in_A = [ index_B_in_A indexB ];
    end
  end
end % find_similar
