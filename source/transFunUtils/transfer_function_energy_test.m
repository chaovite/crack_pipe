function is_valid = transfer_function_energy_test(F)
% This is to test that for wide range of a, the magnitude of reflection coefficient and
% transmission coefficient are both less than 1.
is_valid = false;

for a=[0.25,0.5,1,logspace(-4,4,1000)]
    R = abs(-a*F./(1 + a*F));
    T = abs(1./(1 + a*F));

    if(max(T) > 1 || max(R) > 1)
      return
    end

end        

is_valid = true;
