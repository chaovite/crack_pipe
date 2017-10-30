function is_valid = transfer_function_stability_test(f,F,fmin,fmax)

% Ensure that the phase of the transfer function is bounded by
% + / - pi/2. This means the real part of the transfer function must be
% greater than zero.

[tmp m] = min(abs(f-fmin));
[tmp p] = min(abs(f-fmax));

F = F(m:p);

phase = angle(F);

is_valid = all(abs(phase) <= pi/2);
 
end

