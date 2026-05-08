function v_rot = rotz2d(v, ang)
% Rotate a 3x1 vector in the xy-plane by angle ang about +z.

R = [cos(ang) -sin(ang) 0;
     sin(ang)  cos(ang) 0;
     0         0        1];

v_rot = R * v;
end