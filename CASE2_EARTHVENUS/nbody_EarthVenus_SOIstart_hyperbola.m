function runData = nbody_EarthVenus_SOIstart_hyperbola()
clear; clc; close all;
%  RESTRICTED PRESCRIBED BODY N-BODY: EARTH VENUS
%  COMPARISON VERSION USING ACTUAL DEPARTURE HYPERBOLA STATE

% Common initial state is taken at Earth's SOI exit from the Earth-centred
% departure hyperbola implied by the patched-conic EarthVenus transfer.



% Constants
muS = 1.32712440018e11;   % Sun GM [km^3/s^2]
muE = 3.986004418e5;      % Earth GM [km^3/s^2]
muV = 3.24859e5;          % Venus GM [km^3/s^2]

AU  = 1.495978707e8;      % km
Re  = 6378.1363;          % km
Rv  = 6051.8;             % km

R_E = 1.0*AU;             % Earth heliocentric orbit radius [km]
R_V = 0.723332*AU;        % Venus heliocentric orbit radius [km]

T_E = 365.25*86400;       % s
T_V = 224.701*86400;      % s

nE = 2*pi/T_E;            % rad/s
nV = 2*pi/T_V;            % rad/s

% Hohmann transfer geometry
H = hohmann_transfer(R_E, R_V, muS);

thetaE0 = 0.0;
thetaV0 = mod(thetaE0 + pi - nV*H.tT, 2*pi);

% Earth heliocentric state at departure epoch
rE0 = circular_orbit_pos(R_E, nE, 0, thetaE0);
vE0 = circular_orbit_vel(R_E, nE, 0, thetaE0);

% Heliocentric transfer departure velocity
tHatE = vE0 / norm(vE0);
vD_vec = H.vpT * tHatE;

% Hyperbolic excess velocity relative to Earth
vInf_dep_vec = vD_vec - vE0;
vInf_dep     = norm(vInf_dep_vec);
uInf_dep     = vInf_dep_vec / vInf_dep;

% Earth SOI radius
rSOI_E = R_E * (muE/muS)^(2/5);

% Chosen periapsis radius of Earth departure hyperbola
rp_dep = Re + 300;   % 300 km parking orbit altitude

% ACTUAL Earth-centred hyperbola state at SOI exit
[r_rel_SOI, v_rel_SOI, hyp] = departure_hyperbola_state_at_soi( ...
    muE, rp_dep, vInf_dep_vec, rSOI_E);

% Common heliocentric initial state at Earth SOI exit
r0 = rE0 + r_rel_SOI;
v0 = vE0 + v_rel_SOI;

Y0 = [r0; v0];

fprintf('====================================================\n');
fprintf('PURE RESTRICTED N-BODY: EARTH -> VENUS (SOI START)\n');
fprintf('====================================================\n');
fprintf('Transfer time tT                    = %.6f days\n', H.tT/86400);
fprintf('Initial Earth phase thetaE0         = %.6f rad\n', thetaE0);
fprintf('Initial Venus phase thetaV0         = %.6f rad\n', thetaV0);
fprintf('Earth SOI radius                    = %.6f km\n', rSOI_E);
fprintf('v_inf,dep magnitude                 = %.6f km/s\n', vInf_dep);
fprintf('Departure periapsis rp_dep          = %.6f km\n', rp_dep);
fprintf('Hyperbola eccentricity e_dep        = %.6f\n', hyp.e);
fprintf('Hyperbola beta angle                = %.6f deg\n', rad2deg(hyp.beta));
fprintf('SOI true anomaly theta_SOI          = %.6f deg\n', rad2deg(hyp.theta_soi));
fprintf('Initial heliocentric speed          = %.6f km/s\n', norm(v0));
fprintf('Initial distance from Earth         = %.6f km\n', norm(r0 - rE0));
fprintf('Initial Earth-relative speed        = %.6f km/s\n\n', norm(v_rel_SOI));

%% Parameters for restricted n-body
params.mu       = [muS, muE, muV];
params.R_E      = R_E;
params.R_V      = R_V;
params.nE       = nE;
params.nV       = nV;
params.thetaE0  = thetaE0;
params.thetaV0  = thetaV0;
params.rBodyFun = @rBodies_sun_earth_venusCase3;

% Propagation
dt = 60;                           % s DELTA T
tspan = [0, H.tT + 25*86400];      % enough buffer beyond nominal arrival

tic;
outNB = leapfrog_propagate(@accel_n_body, tspan, Y0, dt, params);
runtimeNB = toc;

tNB = outNB.t;
rNB = outNB.Y(:,1:3)';

% Planet histories
N = numel(tNB);
rEarthHist = zeros(3,N);
rVenusHist = zeros(3,N);

for k = 1:N
    rBodies = rBodies_sun_earth_venusCase3(tNB(k), params);
    rEarthHist(:,k) = rBodies(:,2);
    rVenusHist(:,k) = rBodies(:,3);
end

% Diagnostics
distToEarth = vecnorm(rNB - rEarthHist, 2, 1);
distToVenus = vecnorm(rNB - rVenusHist, 2, 1);

[minDistVenus, idxMin] = min(distToVenus);

fprintf('Runtime                            = %.6f s\n', runtimeNB);
fprintf('Closest approach to Venus          = %.6f km\n', minDistVenus);
fprintf('Time of closest approach           = %.6f days\n', tNB(idxMin)/86400);
fprintf('Nominal Hohmann arrival time       = %.6f days\n', H.tT/86400);
fprintf('Final distance to Venus            = %.6f km\n', distToVenus(end));
fprintf('Final distance to Earth            = %.6f km\n\n', distToEarth(end));

%% Pack outputs for analysis script
runData.out = outNB;
runData.params = params;
runData.H = H;
runData.runtime = runtimeNB;
runData.caseName = 'Earth-Venus restricted n-body (SOI hyperbola start)';
runData.constants.Rv = Rv;

end