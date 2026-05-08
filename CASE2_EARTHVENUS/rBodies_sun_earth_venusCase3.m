function rBodies = rBodies_sun_earth_venusCase3(t, params)
% Sun at origin; Earth and Venus on prescribed circular heliocentric orbits.

rS = [0;0;0];
rE = circular_orbit_pos(params.R_E, params.nE, t, params.thetaE0);
rV = circular_orbit_pos(params.R_V, params.nV, t, params.thetaV0);

rBodies = [rS, rE, rV];
end


