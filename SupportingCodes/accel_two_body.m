function a = accel_two_body(~, r, ~, params)
% accel_two_body  Central gravity acceleration
% a = -mu * r / |r|^3

mu = params.mu;  % km^3/s^2
rn = norm(r);

a = -mu * r / rn^3;
end
