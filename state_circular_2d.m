function [r, v] = state_circular_2d(rmag, mu, theta)
% state_circular_2d  2D circular orbit state at angle theta (rad)
% r = [x;y;0], v tangential (prograde)

r = [rmag*cos(theta); rmag*sin(theta); 0];
vc = sqrt(mu/rmag);

% Tangential unit vector is [-sin, cos]
v = vc * [-sin(theta); cos(theta); 0];
end
