function v = circular_orbit_vel(R, n, t, theta0)
theta = theta0 + n*t;
v = [-R*n*sin(theta); R*n*cos(theta); 0];
end