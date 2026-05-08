function rBodies = rBodies_earth_moonCase1(t, params)
% Earth at origin, Moon on circular orbit around Earth (simplified).
rE = [0;0;0];
rM = circular_orbit_pos(params.R_M, params.n_M, t, params.thetaM0);
rBodies = [rE, rM];
end