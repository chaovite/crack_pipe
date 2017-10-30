function is_valid = transfer_function_dc_test(f,F)
  % The DC-component of the transfer function must be real and positive
  tol = 1e-12;
  is_valid = abs(imag(F(1))) < tol;
  is_valid = and(is_valid,real(F(1))>0);
end
