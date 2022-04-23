#Secuential Jacobi's method implementation
#
function [x_k, k, err] = sec_jacobi(A, b, x_0, iterMax, tol)
  m = rows(A);
  function x_i = get_x_i(i, x_k)
    sum_i = 0;
    for j = 1 : m
      if j != i
        sum_i += A(i,j) * x_k(j);
      end
    end
    x_i = 1/A(i,i) * (b(i) - sum_i);
  endfunction

  x_k = x_0;
  for k = 1 : iterMax
    x_kp1 = zeros(m,1); #initialize x_(k+1)
    for i = 1 : m
      x_kp1(i) = get_x_i(i, x_k);
    end
    x_k = x_kp1;
    err = norm(A * x_k - b);
    if err < tol
      break;
    end
  end
endfunction
