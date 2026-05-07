function a = accel_n_body(t, r, ~, params)
% Acceleration on spacecraft due to multiple bodies with prescribed motion.
% params.mu       : 1xN gravitational parameters (km^3/s^2)
% params.rBodyFun : function handle -> 3xN body positions at time t

rBodies = params.rBodyFun(t, params);   % 3xN
N = numel(params.mu);

a = zeros(3,1);
for i = 1:N
    dr = r - rBodies(:,i);
    dist = norm(dr);
    a = a - params.mu(i) * dr / dist^3;
end
end