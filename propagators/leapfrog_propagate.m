function out = leapfrog_propagate(accelFun, tspan, Y0, dt, params)
% leapfrog_propagate  Velocity-Verlet / Leapfrog integrator for second-order dynamics
% State Y = [r; v] with r,v as 3x1 vectors
% accelFun signature: a = accelFun(t, r, v, params)
% Inputs:
%   accelFun : function handle
%   tspan    : [t0 tf]
%   Y0       : 6x1 initial state
%   dt       : fixed timestep (s)
%   params   : struct (e.g. mu)
% Output struct fields:
%   t : Nx1 time vector
%   Y : Nx6 state history
% Notes:
% - This is a fixed-step symplectic method (good for orbits).
% - Uses kick-drift-kick (velocity Verlet).

t0 = tspan(1);
tf = tspan(2);

% Number of steps obtained by truncation so that the integrator does not step past tf
N = floor((tf - t0)/dt) + 1;
t = (t0 + (0:N-1)'*dt);

Y = zeros(N, 6);
Y(1,:) = Y0(:)';

r = Y0(1:3);
v = Y0(4:6);

% Initial acceleration
a = accelFun(t0, r, v, params);

for k = 1:N-1
    % Kick: half step velocity
    v_half = v + 0.5*dt*a;

    % Drift: full step position
    r_new = r + dt*v_half;

    % New acceleration at updated position
    a_new = accelFun(t(k+1), r_new, v_half, params);

    % Kick: second half step velocity
    v_new = v_half + 0.5*dt*a_new;

    % Store
    Y(k+1,1:3) = r_new';
    Y(k+1,4:6) = v_new';

    % Advance
    r = r_new;
    v = v_new;
    a = a_new;
end

out.t = t;
out.Y = Y;
out.dt = dt;
out.params = params;
end
