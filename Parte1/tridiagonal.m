#Generates a tridiagonal mxm matrix, based on the requirements specifed
#parameters:
#   @p = vector of length m
#   @p = vector of length m
#   @m = size of the square matrix
#output:
#   resultant tridiagonal matrix
function [A] = tridiagonal(p, q, m)
  #Checks if size of vector is m-1
  function is_correct = check_size(v, m)
    is_correct = false;
    if length(v) == m-1
      is_correct = true;
    end
  endfunction


  if m >= 3
    if and(check_size(p, m), check_size(q, m)) #check size of input vectors
      p = [0 p]

      #create diagonal
      d = [2*q(1)];
      for i = 2 : m-1
        d = [d 2*(p(i) + q(i))];
      end
      d = [d 2*p(m)];

      D = diag(d) #diagonal
      Dp1 = diag(q, 1) #first superdiagonal
      p = p(2:m);
      Dm1 = diag(p, -1) #first subdiagonal

      A = Dm1 + D + Dp1;

    else
      error('p and q vector must match the size: m-1')
    end

  else
    error('m must be greater or equal than 3')
  end
endfunction


