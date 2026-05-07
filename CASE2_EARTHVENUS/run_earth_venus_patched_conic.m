function runData = run_earth_venus_patched_conic()
clear; clc; close all;

% running simulation: EARTH-VENUS PATCHED CONIC, with a focus on transfer
% arc
% Same common Earth SOI-exit state as restricted n-body
% Heliocentric cruise propagated with Sun-only gravity

% Constants
muS = 1.32712440018e11;   % km^3/s^2
muE = 3.986004418e5;      % km^3/s^2
muV = 3.24859e5;          % km^3/s^2

AU  = 1.495978707e8;      % km
Re  = 6378.1363;          % km
Rv  = 6051.8;             % km

R_E = 1.0*AU;
R_V = 0.723332*AU;

T_E = 365.25*86400;
T_V = 224.701*86400;

nE = 2*pi/T_E;
nV = 2*pi/T_V;

% Nominal Hohmann geometry
H = hohmann_transfer(R_E, R_V, muS);

thetaE0 = 0.0;
thetaV0 = mod(thetaE0 + pi - nV*H.tT, 2*pi);

% Earth heliocentric departure state
rE0 = circular_orbit_pos(R_E, nE, 0, thetaE0);
vE0 = circular_orbit_vel(R_E, nE, 0, thetaE0);

% Heliocentric transfer departure velocity
tHatE = vE0 / norm(vE0);
vD_vec = H.vpT * tHatE;

% Hyperbolic excess relative to Earth
vInf_dep_vec = vD_vec - vE0;
vInf_dep     = norm(vInf_dep_vec);

% Earth SOI
rSOI_E = R_E * (muE/muS)^(2/5);

% Departure hyperbola choice
rp_dep = Re + 300;   % 300 km parking orbit altitude

[r_rel_SOI, v_rel_SOI, hypDep] = departure_hyperbola_state_at_soi( ...
    muE, rp_dep, vInf_dep_vec, rSOI_E);

% Common heliocentric SOI-exit state
r0 = rE0 + r_rel_SOI;
v0 = vE0 + v_rel_SOI;
Y0 = [r0; v0];
% Venus arrival hyperbola definition 
rSOI_V = R_V * (muV/muS)^(2/5);
rp_arr = Rv + 300;            % 300 km altitude, our design choice

vV_nom = circular_orbit_vel(R_V, nV, H.tT, thetaV0);
tHatV  = vV_nom / norm(vV_nom);
vA_vec = H.vaT * tHatV;       % nominal heliocentric arrival velocity
vInf_arr_vec = vA_vec - vV_nom;
vInf_arr = norm(vInf_arr_vec);

e_arr  = 1 + rp_arr*vInf_arr^2/muV;
vp_arr = sqrt(vInf_arr^2 + 2*muV/rp_arr);

% we Propagate Sun-only heliocentric cruise
dt = 60;
tspan = [0, H.tT + 25*86400];

fprintf('====================================================\n');
fprintf('RUNNING EARTH-VENUS PATCHED-CONIC CRUISE\n');
fprintf('====================================================\n');

tic;
outPC = leapfrog_propagate(@accel_two_body, tspan, Y0, dt, struct('mu',muS));
runtime = toc;

fprintf('Propagation complete. Runtime = %.6f s\n', runtime);

%% Pack outputs
runData.caseName = 'earth_venus_patched_conic';
runData.runtime  = runtime;
runData.out      = outPC;
runData.Y0       = Y0;
runData.dt       = dt;
runData.tspan    = tspan;

runData.constants = struct( ...
    'muS', muS, ...
    'muE', muE, ...
    'muV', muV, ...
    'AU',  AU, ...
    'Re',  Re, ...
    'Rv',  Rv, ...
    'R_E', R_E, ...
    'R_V', R_V, ...
    'nE',  nE, ...
    'nV',  nV);

runData.H        = H;
runData.thetaE0  = thetaE0;
runData.thetaV0  = thetaV0;
runData.rSOI_E   = rSOI_E;
runData.rSOI_V   = rSOI_V;
runData.rp_dep   = rp_dep;
runData.rp_arr   = rp_arr;
runData.hypDep   = hypDep;
runData.vInf_dep = vInf_dep;
runData.vInf_arr = vInf_arr;
runData.e_arr    = e_arr;
runData.vp_arr   = vp_arr;

end