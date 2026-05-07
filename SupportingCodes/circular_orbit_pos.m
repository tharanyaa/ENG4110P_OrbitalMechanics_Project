function r = circular_orbit_pos(R, n, t, theta0)
theta = theta0 + n*t;
r = [R*cos(theta); R*sin(theta); 0];
end