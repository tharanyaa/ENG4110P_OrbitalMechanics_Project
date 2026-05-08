clear; clc; close all;


%CASE STUDY 1: EARTH-MOON : PATCHED CONIC VS RESTRICTED N-BODY
%comparison script for dissertation metrics
% Run both cases
PC = patchedconic_Case1();
NB = NBody_Case1();


% Basic checks

if ~PC.enteredSOI
    error('Patched-conic case did not enter the Moon SOI, so comparison cannot proceed.');
end


% Constants

RM    = PC.constants.RM;
R_M   = PC.constants.R_M;
nM    = PC.constants.nM;
rSOI  = PC.rSOI_M_km;


% PATCHED-CONIC DATA
tE_pc     = PC.tEarth_s(:);
YE_pc     = PC.YEarth;
rSC_E_pc  = YE_pc(:,1:3);
vSC_E_pc  = YE_pc(:,4:6);

rMoon_E_pc = PC.rMoonEarth_km;
vMoon_E_pc = PC.vMoonEarth_kms;
dMoon_E_pc = PC.dMoonEarth_km;

idxSOI_pc  = PC.idxSOI; % Convenience variables
tSOI_pc    = PC.tSOI_s;

tM_pc      = PC.tMoon_s(:);
YM_pc      = PC.YMoon;
rSC_M_pc   = YM_pc(:,1:3);
vSC_M_pc   = YM_pc(:,4:6);
rMagM_pc   = PC.rMagMoon_km;
vMagM_pc   = PC.vMagMoon_kms;

idxCA_pc   = PC.idxClosestApproach;

pcClosestDist = PC.minDistMoon_km;
pcClosestAlt  = PC.minAltMoon_km;
pcClosestTime = PC.tSOI_s + tM_pc(idxCA_pc);
pcClosestVel  = vMagM_pc(idxCA_pc);

runtimePC = PC.runtime;
memoryPC_MB = PC.memory_MB;

% Fair step counts for actual patched trajectory
stepsPC_earth = idxSOI_pc;
stepsPC_moon  = numel(tM_pc);
stepsPC_total = stepsPC_earth + stepsPC_moon - 1;

% RESTRICTED N-BODY DATA
outNB   = NB.out;
params  = NB.params;

tNB = outNB.t(:);
YNB = outNB.Y;

rSC_NB = YNB(:,1:3);
vSC_NB = YNB(:,4:6);

N = numel(tNB);

rMoon_NB = zeros(N,3);
vMoon_NB = zeros(N,3);

for k = 1:N
    rBodies = params.rBodyFun(tNB(k), params);
    rMoon_NB(k,:) = rBodies(:,2).';

    thetaM = params.thetaM0 + params.n_M*tNB(k);
    vMoon_NB(k,:) = [-R_M*nM*sin(thetaM), ...
                      R_M*nM*cos(thetaM), ...
                      0];
end

rRelMoon_NB = rSC_NB - rMoon_NB;
vRelMoon_NB = vSC_NB - vMoon_NB;

distMoon_NB = vecnorm(rRelMoon_NB, 2, 2);
relSpeedMoon_NB = vecnorm(vRelMoon_NB, 2, 2);

[minDistMoon_NB, idxMinNB] = min(distMoon_NB);
minAltMoon_NB = minDistMoon_NB - RM;

nbClosestDist = minDistMoon_NB;
nbClosestAlt  = minAltMoon_NB;
nbClosestTime = tNB(idxMinNB);
nbClosestVel  = relSpeedMoon_NB(idxMinNB);

runtimeNB = NB.runtime;
memoryNB_MB = NB.memory_MB;
stepsNB   = numel(tNB);

% First n-body crossing of same SOI threshold
idxSOI_nb = find(distMoon_NB <= rSOI, 1, 'first');

if isempty(idxSOI_nb)
    tSOI_nb = NaN;
    soiTimingDiff_days = NaN;
else
    tSOI_nb = tNB(idxSOI_nb);
    soiTimingDiff_days = (tSOI_nb - tSOI_pc)/86400;
end

% Supporting diagnostic only
nbFinalSep = distMoon_NB(end);

% COMMON-INTERVAL TRAJECTORY COMPARISONS

% 1. Earth-phase inertial position and velocity deviation
% Compare over common Earth-phase interval:
% from launch to patched-conic SOI entry
tEndEarth = min(tSOI_pc, tNB(end));

maskNB_E = tNB <= tEndEarth;
maskPC_E = tE_pc <= tEndEarth;

tCommonEarth = tNB(maskNB_E);

rPC_E_interp = interp1(tE_pc(maskPC_E), rSC_E_pc(maskPC_E,:), tCommonEarth, 'linear');
vPC_E_interp = interp1(tE_pc(maskPC_E), vSC_E_pc(maskPC_E,:), tCommonEarth, 'linear');

deltaR_earth = vecnorm(rSC_NB(maskNB_E,:) - rPC_E_interp, 2, 2);
deltaV_earth = vecnorm(vSC_NB(maskNB_E,:) - vPC_E_interp, 2, 2);

earthPos_max_km   = max(deltaR_earth);
earthPos_rms_km   = sqrt(mean(deltaR_earth.^2));
earthPos_final_km = deltaR_earth(end);

earthVel_max_kmps   = max(deltaV_earth);
earthVel_rms_kmps   = sqrt(mean(deltaV_earth.^2));
earthVel_final_kmps = deltaV_earth(end);

% 2. Moon-relative position deviation over common approach interval
% Compare Moon-relative geometry from launch to end of patched-conic run
tEndCommon = min(tNB(end), tSOI_pc + tM_pc(end));

maskNB_C = tNB <= tEndCommon;

% Earth phase: only keep launch -> SOI entry
tE_pc_toSOI       = tE_pc(1:idxSOI_pc);
rRelMoon_PC_earth = rSC_E_pc(1:idxSOI_pc,:) - rMoon_E_pc(1:idxSOI_pc,:);
vRelMoon_PC_earth = vSC_E_pc(1:idxSOI_pc,:) - vMoon_E_pc(1:idxSOI_pc,:);

% Moon phase already Moon-relative
rRelMoon_PC_moon = rSC_M_pc;
vRelMoon_PC_moon = vSC_M_pc;

% Continuous stitched patched-conic Moon-relative history
tPC_full = [tE_pc_toSOI; tSOI_pc + tM_pc(2:end)];
rRelMoon_PC_full = [rRelMoon_PC_earth; rRelMoon_PC_moon(2:end,:)];
vRelMoon_PC_full = [vRelMoon_PC_earth; vRelMoon_PC_moon(2:end,:)];

% Safety guard against duplicate times
[tPC_full, ia] = unique(tPC_full, 'stable');
rRelMoon_PC_full = rRelMoon_PC_full(ia,:);
vRelMoon_PC_full = vRelMoon_PC_full(ia,:);

tCommonMoon = tNB(maskNB_C);

rRelMoon_PC_interp = interp1(tPC_full, rRelMoon_PC_full, tCommonMoon, 'linear');
vRelMoon_PC_interp = interp1(tPC_full, vRelMoon_PC_full, tCommonMoon, 'linear');

deltaR_relMoon = vecnorm(rRelMoon_NB(maskNB_C,:) - rRelMoon_PC_interp, 2, 2);
deltaV_relMoon = vecnorm(vRelMoon_NB(maskNB_C,:) - vRelMoon_PC_interp, 2, 2);

moonRelPos_max_km   = max(deltaR_relMoon);
moonRelPos_rms_km   = sqrt(mean(deltaR_relMoon.^2));
moonRelPos_final_km = deltaR_relMoon(end);

moonRelVel_max_kmps   = max(deltaV_relMoon);
moonRelVel_rms_kmps   = sqrt(mean(deltaV_relMoon.^2));
moonRelVel_final_kmps = deltaV_relMoon(end);

% 3. Stitched actual-conic Moon-distance history
% This is the ACTUAL patched-conic mission distance history
tPC_dist = [tE_pc_toSOI; tSOI_pc + tM_pc(2:end)];
dMoon_PC_stitched = [dMoon_E_pc(1:idxSOI_pc); rMagM_pc(2:end)];

[tPC_dist, ia_dist] = unique(tPC_dist, 'stable');
dMoon_PC_stitched = dMoon_PC_stitched(ia_dist);

% ENCOUNTER METRICS
timingShift_days = (nbClosestTime - pcClosestTime)/86400;

distanceError_km = nbClosestDist - pcClosestDist;
altitudeError_km = nbClosestAlt  - pcClosestAlt;
speedError_kmps  = nbClosestVel  - pcClosestVel;

% PRINT COMPARISON SUMMARY
fprintf('====================================================\n');
fprintf('REVISED COMPARISON: EARTH-MOON PATCHED CONIC VS RESTRICTED N-BODY\n');
fprintf('====================================================\n\n');

fprintf('[SOI Entry Comparison]\n');
fprintf('Patched-conic SOI entry time            = %.6f days since launch\n', tSOI_pc/86400);
if isnan(tSOI_nb)
    fprintf('Restricted n-body SOI threshold crossing = not reached\n');
    fprintf('SOI entry timing difference             = n/a\n\n');
else
    fprintf('Restricted n-body SOI threshold crossing = %.6f days since launch\n', tSOI_nb/86400);
    fprintf('SOI entry timing difference             = %.6f days\n\n', soiTimingDiff_days);
end

fprintf('[Closest Approach Timing]\n');
fprintf('Patched-conic closest approach time     = %.6f days since launch\n', pcClosestTime/86400);
fprintf('Restricted n-body closest approach time = %.6f days since launch\n', nbClosestTime/86400);
fprintf('Timing shift                            = %.6f days\n\n', timingShift_days);

fprintf('[Closest Approach Geometry]\n');
fprintf('Patched-conic minimum Moon distance     = %.6f km\n', pcClosestDist);
fprintf('Restricted n-body minimum Moon distance = %.6f km\n', nbClosestDist);
fprintf('Distance difference                     = %.6f km\n\n', distanceError_km);

fprintf('Patched-conic minimum altitude          = %.6f km\n', pcClosestAlt);
fprintf('Restricted n-body minimum altitude      = %.6f km\n', nbClosestAlt);
fprintf('Altitude difference                     = %.6f km\n\n', altitudeError_km);

fprintf('[Relative Speed at Closest Approach]\n');
fprintf('Patched-conic Moon-relative speed       = %.6f km/s\n', pcClosestVel);
fprintf('Restricted n-body Moon-relative speed   = %.6f km/s\n', nbClosestVel);
fprintf('Speed difference                        = %.6f km/s\n\n', speedError_kmps);

fprintf('[Earth-Phase Inertial Deviation]\n');
fprintf('Max position deviation                  = %.6f km\n', earthPos_max_km);
fprintf('RMS position deviation                  = %.6f km\n', earthPos_rms_km);
fprintf('Final position deviation                = %.6f km\n\n', earthPos_final_km);

fprintf('Max velocity deviation                  = %.6f km/s\n', earthVel_max_kmps);
fprintf('RMS velocity deviation                  = %.6f km/s\n', earthVel_rms_kmps);
fprintf('Final velocity deviation                = %.6f km/s\n\n', earthVel_final_kmps);

fprintf('[Moon-Relative Deviation Over Common Run]\n');
fprintf('Max Moon-relative position deviation    = %.6f km\n', moonRelPos_max_km);
fprintf('RMS Moon-relative position deviation    = %.6f km\n', moonRelPos_rms_km);
fprintf('Final Moon-relative position deviation  = %.6f km\n\n', moonRelPos_final_km);

fprintf('Max Moon-relative velocity deviation    = %.6f km/s\n', moonRelVel_max_kmps);
fprintf('RMS Moon-relative velocity deviation    = %.6f km/s\n', moonRelVel_rms_kmps);
fprintf('Final Moon-relative velocity deviation  = %.6f km/s\n\n', moonRelVel_final_kmps);

fprintf('[Computational Cost]\n');
fprintf('Patched-conic runtime                   = %.6f s\n', runtimePC);
fprintf('Restricted n-body runtime               = %.6f s\n', runtimeNB);
fprintf('Patched-conic stored-array memory       = %.6f MB\n', memoryPC_MB);
fprintf('Restricted n-body stored-array memory   = %.6f MB\n', memoryNB_MB);
fprintf('Patched-conic Earth-phase steps         = %d\n', stepsPC_earth);
fprintf('Patched-conic Moon-phase steps          = %d\n', stepsPC_moon);
fprintf('Patched-conic total steps               = %d\n', stepsPC_total);
fprintf('Restricted n-body total steps           = %d\n\n', stepsNB);

fprintf('[Supporting Diagnostic]\n');
fprintf('Restricted n-body final spacecraft-Moon separation = %.6f km\n\n', nbFinalSep);

% PLOTS

% Plot 1: Spacecraft-Moon distance comparison
% Uses the ACTUAL stitched patched-conic trajectory
figure;
plot(tPC_dist/3600, dMoon_PC_stitched, 'LineWidth', 1.4); hold on;
plot(tNB/3600, distMoon_NB, 'LineWidth', 1.4);
yline(rSOI, '--');
grid on;
xlabel('Time since launch (hours)');
ylabel('Spacecraft-Moon distance (km)');
title('Spacecraft-Moon distance: patched-conic vs restricted n-body');
legend('Patched-conic','Restricted n-body','Moon SOI radius','Location','best');

% Plot 2: Earth-phase inertial position deviation
figure;
plot(tCommonEarth/3600, deltaR_earth, 'LineWidth', 1.4);
grid on;
xlabel('Time since launch (hours)');
ylabel('||r_{NB} - r_{PC}|| (km)');
title('Earth-phase inertial position deviation');

% Plot 3: Earth-phase inertial velocity deviation
figure;
plot(tCommonEarth/3600, deltaV_earth, 'LineWidth', 1.4);
grid on;
xlabel('Time since launch (hours)');
ylabel('||v_{NB} - v_{PC}|| (km/s)');
title('Earth-phase inertial velocity deviation');

% Plot 4: Moon-relative position deviation
figure;
plot(tCommonMoon/3600, deltaR_relMoon, 'LineWidth', 1.4);
grid on;
xlabel('Time since launch (hours)');
ylabel('Moon-relative position deviation (km)');
title('Moon-relative position deviation: patched-conic vs restricted n-body');

% Plot 5: Moon-relative velocity deviation
figure;
plot(tCommonMoon/3600, deltaV_relMoon, 'LineWidth', 1.4);
grid on;
xlabel('Time since launch (hours)');
ylabel('Moon-relative velocity deviation (km/s)');
title('Moon-relative velocity deviation: patched-conic vs restricted n-body');

%% Plot 6: Patched-conic Moon-centred phase
figure;
plot(rSC_M_pc(:,1), rSC_M_pc(:,2), 'LineWidth', 1.4); hold on;
th = linspace(0,2*pi,400);
plot(RM*cos(th), RM*sin(th), '--', 'LineWidth', 1.2);
plot(0,0,'k.','MarkerSize',16);
plot(rSC_M_pc(idxCA_pc,1), rSC_M_pc(idxCA_pc,2), 'o', 'MarkerSize', 8, 'LineWidth', 1.4);
axis equal; grid on;
xlabel('x_M (km)');
ylabel('y_M (km)');
title('Patched-conic Moon-centred phase');
legend('Spacecraft','Moon surface','Moon centre','Closest approach','Location','best');

% Plot 7: Restricted n-body Moon-relative encounter geometry
figure;
plot(rRelMoon_NB(:,1), rRelMoon_NB(:,2), 'LineWidth', 1.4); hold on;
th = linspace(0,2*pi,400);
plot(RM*cos(th), RM*sin(th), '--', 'LineWidth', 1.2);
plot(0,0,'k.','MarkerSize',16);
plot(rRelMoon_NB(idxMinNB,1), rRelMoon_NB(idxMinNB,2), 'o', 'MarkerSize', 8, 'LineWidth', 1.4);
axis equal; grid on;
xlabel('x_M (km)');
ylabel('y_M (km)');
title('Restricted n-body Moon-relative encounter geometry');
legend('Spacecraft','Moon surface','Moon centre','Closest approach','Location','best');