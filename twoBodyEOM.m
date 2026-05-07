function dYdt = twoBodyEOM(~, Y, mu)
% twoBodyEOM: 2-body equations of motion for a spacecraft about a central body.
% State Y = [r; v] where r and v are 3x1 vectors.

r = Y(1:3);
v = Y(4:6);

rnorm = norm(r);
a = -mu * r / rnorm^3;

dYdt = [v; a];
end
