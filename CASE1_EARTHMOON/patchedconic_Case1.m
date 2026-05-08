function runData = patchedconic_Case1()
clear; clc; close all;

% PATCHED-CONIC model EARTH-MOON CASE
% Target-matching Hohmann-style translunar initialisation

% Constants
muE = 3.986004418e5;      % Earth gravitational parameter [km^3/s^2]
muM = 4.9048695e3;        % Moon gravitational parameter  [km^3/s^2]
Re  = 6378.1363;          % Earth radius [km]
RM  = 1737.4;             % Moon radius  [km]
R_M = 384400;             % Mean Earth-Moon distance [km]

% Moon mean motion 
nM = 2*pi/(27.321661*86400);   % [rad/s]

% Lunar sphere of influence (SOI)
rSOI_M = R_M * (muM/muE)^(2/5);

% Parking orbit / transfer setup
r1 = Re + 500;                 % parking orbit radius [km]
n1 = sqrt(muE / r1^3);         % parking orbit mean motion [rad/s]

H = hohmann_transfer(r1, R_M, muE);

if ~isfield(H,'vaT')
    H.vaT = sqrt(muE * (2/R_M - 1/H.aT));
end

fprintf('====================================================\n');
fprintf('SIMPLE PATCHED-CONIC TRANSLUNAR INITIALISATION\n');
fprintf('====================================================\n');
fprintf('Transfer semi-major axis aT   = %.6f km\n', H.aT);
fprintf('Transfer time tT              = %.6f days\n', H.tT/86400);
fprintf('Perigee transfer speed vpT    = %.6f km/s\n', H.vpT);
fprintf('Apogee transfer speed vaT     = %.6f km/s\n', H.vaT);
fprintf('Departure circular speed vc1  = %.6f km/s\n', H.vc1);
fprintf('Initial delta-v dv1           = %.6f km/s\n', H.dv1);
fprintf('Moon SOI radius               = %.6f km\n', rSOI_M);

% Target-matching setup
thetaSC0 = 0;

rng(1);
thetaM0_abs = 2*pi*rand();   % Moon phase at absolute time t = 0

deltaAim = -0.2;  % radians  – same value as n-body case
tL = solve_launch_time_hohmann(thetaSC0, thetaM0_abs, n1, nM, H.tT, deltaAim);
tArr = tL + H.tT;

fprintf('Moon initial phase thetaM0    = %.6f rad\n', thetaM0_abs);
fprintf('Solved launch time tL         = %.6f days\n', tL/86400);
fprintf('Arrival time tArr             = %.6f days\n', tArr/86400);

% Moon phase at launch epoch
thetaM0_launch = mod(thetaM0_abs + nM*tL, 2*pi);

% Spacecraft state at launch
thetaSC_launch = thetaSC0 + n1*tL;
[r0E, vCirc0] = state_circular_2d(r1, muE, thetaSC_launch);

vHat = vCirc0 / norm(vCirc0);
v0E = H.vpT * vHat;

Y0E = [r0E; v0E];

% Propagate Earth-centred phase
tic;

dt = 2;
outE = leapfrog_propagate(@accel_two_body, [0, H.tT], Y0E, dt, struct('mu',muE));

tE = outE.t(:);
YE = outE.Y;

rSC_E = YE(:,1:3);
vSC_E = YE(:,4:6);

N = numel(tE);
rMoon_E = zeros(N,3);
vMoon_E = zeros(N,3);
rRelMoon_E = zeros(N,3);
vRelMoon_E = zeros(N,3);
dMoon_E = zeros(N,1);

for k = 1:N
    thetaM = thetaM0_launch + nM*tE(k);

    rMoon_E(k,:) = [R_M*cos(thetaM), R_M*sin(thetaM), 0];
    vMoon_E(k,:) = [-R_M*nM*sin(thetaM), R_M*nM*cos(thetaM), 0];

    rRelMoon_E(k,:) = rSC_E(k,:) - rMoon_E(k,:);
    vRelMoon_E(k,:) = vSC_E(k,:) - vMoon_E(k,:);
    dMoon_E(k) = norm(rRelMoon_E(k,:));
end

% Find SOI entry
idxSOI = find(dMoon_E <= rSOI_M, 1, 'first');

fprintf('\n====================================================\n');
fprintf('PATCHED-CONIC EARTH PHASE\n');
fprintf('====================================================\n');

if isempty(idxSOI)
    fprintf('Moon SOI was NOT entered during the Earth-centred phase.\n');

    [minDistEarthPhase, idxMinE] = min(dMoon_E);
    fprintf('Closest Moon distance during Earth phase = %.6f km\n', minDistEarthPhase);
    fprintf('Closest Moon altitude during Earth phase = %.6f km\n', minDistEarthPhase - RM);

    runtime_s = toc;

    %  MATLAB memory footprint
    tE_mem = tE;
    YE_mem = YE;
    rSC_E_mem = rSC_E;
    vSC_E_mem = vSC_E;
    rMoon_E_mem = rMoon_E;
    vMoon_E_mem = vMoon_E;
    rRelMoon_E_mem = rRelMoon_E;
    vRelMoon_E_mem = vRelMoon_E;
    dMoon_E_mem = dMoon_E;
    Y0E_mem = Y0E;

    varsPC = whos('tE_mem','YE_mem','rSC_E_mem','vSC_E_mem','rMoon_E_mem','vMoon_E_mem', ...
                  'rRelMoon_E_mem','vRelMoon_E_mem','dMoon_E_mem','Y0E_mem');
    memory_bytes = sum([varsPC.bytes]);
    memory_MB = memory_bytes / 1024^2;

    fprintf('Indicative simulation memory footprint = %.6f MB\n', memory_MB);

    runData = struct();
    runData.enteredSOI = false;
    runData.tEarth_s = tE;
    runData.YEarth = YE;
    runData.rMoonEarth_km = rMoon_E;
    runData.dMoonEarth_km = dMoon_E;
    runData.constants.muE = muE;
    runData.constants.muM = muM;
    runData.constants.Re = Re;
    runData.constants.RM = RM;
    runData.constants.R_M = R_M;
    runData.constants.nM = nM;
    runData.rSOI_M_km = rSOI_M;
    runData.H = H;
    runData.runtime = runtime_s;
    runData.memory_bytes = memory_bytes;
    runData.memory_MB = memory_MB;
    return;
end

tSOI = tE(idxSOI);

fprintf('Moon SOI entered.\n');
fprintf('SOI entry time                = %.6f days since launch\n', tSOI/86400);
fprintf('Moon distance at entry        = %.6f km\n', dMoon_E(idxSOI));

% Patch-point Moon-relative state
r0M = rRelMoon_E(idxSOI,:).';
v0M = vRelMoon_E(idxSOI,:).';

fprintf('Moon-relative speed at entry  = %.6f km/s\n', norm(v0M));

% Propagate Moon-centred phase
tRemain = H.tT - tSOI;
outM = leapfrog_propagate(@accel_two_body, [0, tRemain], [r0M; v0M], dt, struct('mu',muM));

tM = outM.t(:);
YM = outM.Y;

rSC_M = YM(:,1:3);
vSC_M = YM(:,4:6);

rMagM = sqrt(sum(rSC_M.^2,2));
vMagM = sqrt(sum(vSC_M.^2,2));

% Check for lunar impact
idxImpact = find(rMagM <= RM, 1, 'first');

fprintf('\n====================================================\n');
fprintf('PATCHED-CONIC MOON PHASE\n');
fprintf('====================================================\n');

if isempty(idxImpact)
    [minDistMoon, idxCA] = min(rMagM);

    fprintf('No lunar impact detected in Moon-centred phase.\n');
    fprintf('Closest approach time         = %.6f days since SOI entry\n', tM(idxCA)/86400);
    fprintf('Minimum Moon distance         = %.6f km\n', minDistMoon);
    fprintf('Minimum Moon altitude         = %.6f km\n', minDistMoon - RM);
    fprintf('Moon-relative speed at CA     = %.6f km/s\n', vMagM(idxCA));

    impactFlag = false;
else
    fprintf('Lunar impact detected.\n');
    fprintf('Impact time after SOI entry   = %.6f days\n', tM(idxImpact)/86400);
    fprintf('Impact distance               = %.6f km\n', rMagM(idxImpact));
    fprintf('Impact speed                  = %.6f km/s\n', vMagM(idxImpact));

    idxCA = idxImpact;
    minDistMoon = rMagM(idxImpact);
    impactFlag = true;
end

% Final diagnostic values
rFinalM = rSC_M(end,:);
vFinalM = vSC_M(end,:);

fprintf('\n[Final Moon-centred State]\n');
fprintf('Final Moon-relative distance  = %.6f km\n', norm(rFinalM));
fprintf('Final Moon-relative speed     = %.6f km/s\n', norm(vFinalM));

runtime_s = toc;

%  MATLAB memory footprint
tE_mem = tE;
YE_mem = YE;
rSC_E_mem = rSC_E;
vSC_E_mem = vSC_E;
rMoon_E_mem = rMoon_E;
vMoon_E_mem = vMoon_E;
rRelMoon_E_mem = rRelMoon_E;
vRelMoon_E_mem = vRelMoon_E;
dMoon_E_mem = dMoon_E;
tM_mem = tM;
YM_mem = YM;
rSC_M_mem = rSC_M;
vSC_M_mem = vSC_M;
rMagM_mem = rMagM;
vMagM_mem = vMagM;
Y0E_mem = Y0E;
r0M_mem = r0M;
v0M_mem = v0M;

varsPC = whos('tE_mem','YE_mem','rSC_E_mem','vSC_E_mem','rMoon_E_mem','vMoon_E_mem', ...
              'rRelMoon_E_mem','vRelMoon_E_mem','dMoon_E_mem', ...
              'tM_mem','YM_mem','rSC_M_mem','vSC_M_mem','rMagM_mem','vMagM_mem', ...
              'Y0E_mem','r0M_mem','v0M_mem');
memory_bytes = sum([varsPC.bytes]);
memory_MB = memory_bytes / 1024^2;

fprintf('Indicative simulation memory footprint = %.6f MB\n', memory_MB);

% Store results in runData for use for comparison later
runData = struct();

runData.enteredSOI = true;
runData.impactFlag = impactFlag;
runData.runtime = runtime_s;
runData.memory_bytes = memory_bytes;
runData.memory_MB = memory_MB;

runData.tEarth_s = tE;
runData.YEarth = YE;
runData.rMoonEarth_km = rMoon_E;
runData.vMoonEarth_kms = vMoon_E;
runData.dMoonEarth_km = dMoon_E;

runData.idxSOI = idxSOI;
runData.tSOI_s = tSOI;
runData.rSOI_M_km = rSOI_M;
runData.rPatchMoon_km = r0M;
runData.vPatchMoon_kms = v0M;

runData.tMoon_s = tM;
runData.YMoon = YM;
runData.rMagMoon_km = rMagM;
runData.vMagMoon_kms = vMagM;

runData.idxClosestApproach = idxCA;
runData.minDistMoon_km = minDistMoon;
runData.minAltMoon_km = minDistMoon - RM;

if impactFlag
    runData.idxImpact = idxImpact;
    runData.tImpact_s = tM(idxImpact);
else
    runData.idxImpact = [];
    runData.tImpact_s = [];
end

runData.constants.muE = muE;
runData.constants.muM = muM;
runData.constants.Re = Re;
runData.constants.RM = RM;
runData.constants.R_M = R_M;
runData.constants.nM = nM;

runData.thetaSC0 = thetaSC0;
runData.thetaM0_abs = thetaM0_abs;
runData.thetaM0_launch = thetaM0_launch;
runData.tL = tL;
runData.tArr = tArr;
runData.thetaSC_launch = thetaSC_launch;
runData.H = H;

end