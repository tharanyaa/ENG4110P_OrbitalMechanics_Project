function runData = NBody_Case1()
clear; clc; close all;

% RESTRICTED EARTH-MOON N-BODY CASE
% Target-matching Hohmann-style translunar initialisation

% Constants
muE = 3.986004418e5;
muM = 4.9048695e3;
Re  = 6378.1363;
RM  = 1737.4;
R_M = 384400;

nM = 2*pi/(27.321661*86400);

% Parking orbit and  transfer setup
r1 = Re + 500;
n1 = sqrt(muE / r1^3);

H = hohmann_transfer(r1, R_M, muE);

if ~isfield(H, 'vaT')
    H.vaT = sqrt(muE * (2/R_M - 1/H.aT));
end

fprintf('====================================================\n');
fprintf('Target-matching translunar Hohmann-style initialisation\n');
fprintf('====================================================\n');
fprintf('Transfer semi-major axis aT   = %.6f km\n', H.aT);
fprintf('Transfer time tT              = %.6f days\n', H.tT/86400);
fprintf('Perigee transfer speed vpT    = %.6f km/s\n', H.vpT);
fprintf('Apogee transfer speed vaT     = %.6f km/s\n', H.vaT);
fprintf('Departure circular speed vc1  = %.6f km/s\n', H.vc1);
fprintf('Initial delta-v dv1           = %.6f km/s\n', H.dv1);

% -----------------------------
% Target-matching setup
% -----------------------------
thetaSC0 = 0;

rng(1);
thetaM0 = 2*pi*rand();


deltaAim = -0.2;  % radians (~0.3 deg) – adjust if needed
tL = solve_launch_time_hohmann(thetaSC0, thetaM0, n1, nM, H.tT, deltaAim);
tArr = tL + H.tT;

fprintf('Moon initial phase thetaM0    = %.6f rad\n', thetaM0);
fprintf('Solved launch time tL         = %.6f days\n', tL/86400);
fprintf('Arrival time tArr             = %.6f days\n', tArr/86400);

% -----------------------------
% Spacecraft state at launch
% -----------------------------
thetaSC_launch = thetaSC0 + n1*tL;

[rL, vCircL] = state_circular_2d(r1, muE, thetaSC_launch);

vHat = vCircL / norm(vCircL);
vL_transfer = H.vpT * vHat;

Y0 = [rL; vL_transfer];

% Simulation settings
tspan = [0, H.tT];
dt = 2;

% Restricted n-body parameters
params.mu = [muE muM];
params.R_M = R_M;
params.n_M = nM;
params.thetaM0 = mod(thetaM0 + nM*tL, 2*pi);
params.rBodyFun = @rBodies_earth_moonCase1;

% Propagation
fprintf('\n====================================================\n');
fprintf('Running restricted n-body propagation from launch\n');
fprintf('====================================================\n');

tic;
out = leapfrog_propagate(@accel_n_body, tspan, Y0, dt, params);
runtime_s = toc;

fprintf('Propagation complete.\n');
fprintf('Runtime = %.6f s\n', runtime_s);

% Simple phase and geometry checks
thetaSC_arr_ideal = thetaSC_launch + pi;
thetaM_arr_ideal  = thetaM0 + nM*tArr;
phase_err = wrapToPi(thetaSC_arr_ideal - thetaM_arr_ideal);

fprintf('\n====================================================\n');
fprintf('Target-matching checks\n');
fprintf('====================================================\n');
fprintf('Ideal phase error at arrival = %.6e rad\n', phase_err);

rSC_end = out.Y(end,1:3).';
rBodies_end = params.rBodyFun(H.tT, params);
rMoon_end = rBodies_end(:,2);

sep_end = norm(rSC_end - rMoon_end);

fprintf('Spacecraft-Moon separation at arrival = %.6f km\n', sep_end);

%  MATLAB memory footprint
tNB_mem = out.t(:);
YNB_mem = out.Y;
rSC_end_mem = rSC_end;
rBodies_end_mem = rBodies_end;
rMoon_end_mem = rMoon_end;
Y0_mem = Y0;

varsNB = whos('tNB_mem','YNB_mem','rSC_end_mem','rBodies_end_mem','rMoon_end_mem','Y0_mem');
memory_bytes = sum([varsNB.bytes]);
memory_MB = memory_bytes / 1024^2;

fprintf('Indicative simulation memory footprint = %.6f MB\n', memory_MB);

% Packing outputs for analysis script
runData.caseName = 'earth_moon_restricted_nbody';
runData.out = out;
runData.params = params;
runData.H = H;
runData.runtime = runtime_s;
runData.memory_bytes = memory_bytes;
runData.memory_MB = memory_MB;

runData.constants.muE = muE;
runData.constants.muM = muM;
runData.constants.Re  = Re;
runData.constants.RM  = RM;
runData.constants.R_M = R_M;
runData.constants.nM  = nM;

runData.r1 = r1;
runData.n1 = n1;
runData.thetaSC0 = thetaSC0;
runData.thetaM0 = thetaM0;
runData.thetaSC_launch = thetaSC_launch;
runData.tL = tL;
runData.tArr = tArr;
runData.phase_err = phase_err;
runData.sep_end = sep_end;
runData.Y0 = Y0;
runData.dt = dt;
runData.tspan = tspan;

end