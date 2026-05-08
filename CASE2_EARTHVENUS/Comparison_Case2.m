clear; clc; close all;

% Run both models
R_nb = nbody_EarthVenus_SOIstart_hyperbola();
R_pc = run_earth_venus_patched_conic();

% EXTRACT N-BODY OUTPUTS
outNB    = R_nb.out;
paramsNB = R_nb.params;
H_nb     = R_nb.H;
Rv_nb    = R_nb.constants.Rv;

tNB = outNB.t(:);
YNB = outNB.Y;
rNB = YNB(:,1:3);
vNB = YNB(:,4:6);

Nnb = numel(tNB);

rEarthNB = zeros(Nnb,3);
rVenusNB = zeros(Nnb,3);
vVenusNB = zeros(Nnb,3);

for k = 1:Nnb
    rBodies = paramsNB.rBodyFun(tNB(k), paramsNB);
    rEarthNB(k,:) = rBodies(:,2)';
    rVenusNB(k,:) = rBodies(:,3)';
    vVenusNB(k,:) = circular_orbit_vel(paramsNB.R_V, paramsNB.nV, tNB(k), paramsNB.thetaV0)';
end

distEarthNB = vecnorm(rNB - rEarthNB, 2, 2);
distVenusNB = vecnorm(rNB - rVenusNB, 2, 2);
relVelVenusNB = vecnorm(vNB - vVenusNB, 2, 2);

[minDistVenusNB, idxMinNB] = min(distVenusNB);
minAltVenusNB = minDistVenusNB - Rv_nb;

A_nb.runtime = R_nb.runtime;
A_nb.closestApproach.timeDays     = tNB(idxMinNB)/86400;
A_nb.closestApproach.distanceKm   = minDistVenusNB;
A_nb.closestApproach.altitudeKm   = minAltVenusNB;
A_nb.closestApproach.relSpeedKmps = relVelVenusNB(idxMinNB);
A_nb.closestApproach.idx          = idxMinNB;
A_nb.nominalHohmannArrivalDays    = H_nb.tT/86400;

% EXTRACT PATCHED-CONIC OUTPUTS
=outPC   = R_pc.out;
H_pc    = R_pc.H;
Rv_pc   = R_pc.constants.Rv;
R_E_pc  = R_pc.constants.R_E;
R_V_pc  = R_pc.constants.R_V;
nE_pc   = R_pc.constants.nE;
nV_pc   = R_pc.constants.nV;
thetaE0 = R_pc.thetaE0;
thetaV0 = R_pc.thetaV0;

tPC = outPC.t(:);
YPC = outPC.Y;
rPC = YPC(:,1:3);
vPC = YPC(:,4:6);

Npc = numel(tPC);

rEarthPC = zeros(Npc,3);
rVenusPC = zeros(Npc,3);
vVenusPC = zeros(Npc,3);

for k = 1:Npc
    rEarthPC(k,:) = circular_orbit_pos(R_E_pc, nE_pc, tPC(k), thetaE0)';
    rVenusPC(k,:) = circular_orbit_pos(R_V_pc, nV_pc, tPC(k), thetaV0)';
    vVenusPC(k,:) = circular_orbit_vel(R_V_pc, nV_pc, tPC(k), thetaV0)';
end

distEarthPC = vecnorm(rPC - rEarthPC, 2, 2);
distVenusPC = vecnorm(rPC - rVenusPC, 2, 2);
relVelVenusPC = vecnorm(vPC - vVenusPC, 2, 2);

[minDistVenusPC, idxMinPC] = min(distVenusPC);
minAltVenusPC = minDistVenusPC - Rv_pc;

A_pc.runtime = R_pc.runtime;
A_pc.closestApproach.timeDays     = tPC(idxMinPC)/86400;
A_pc.closestApproach.distanceKm   = minDistVenusPC;
A_pc.closestApproach.altitudeKm   = minAltVenusPC;
A_pc.closestApproach.relSpeedKmps = relVelVenusPC(idxMinPC);
A_pc.closestApproach.idx          = idxMinPC;

A_pc.nominalArrival.periapsisRadiusKm   = R_pc.rp_arr;
A_pc.nominalArrival.periapsisAltitudeKm = R_pc.rp_arr - Rv_pc;
A_pc.nominalArrival.vInfKmps            = R_pc.vInf_arr;
A_pc.nominalArrival.periapsisSpeedKmps  = R_pc.vp_arr;
A_pc.nominalArrival.hohmannArrivalDays  = H_pc.tT/86400;

% DIRECT PATCHED vs N-BODY CLOSEST-APPROACH COMPARISON
timeDiffClosest_days     = A_nb.closestApproach.timeDays     - A_pc.closestApproach.timeDays;
distanceDiffClosest_km   = A_nb.closestApproach.distanceKm   - A_pc.closestApproach.distanceKm;
altitudeDiffClosest_km   = A_nb.closestApproach.altitudeKm   - A_pc.closestApproach.altitudeKm;
relSpeedDiffClosest_kmps = A_nb.closestApproach.relSpeedKmps - A_pc.closestApproach.relSpeedKmps;

% NOMINAL PATCHED ARRIVAL TARGET vs N-BODY REALISED OUTCOME
timingShiftVsNominal_days = A_nb.closestApproach.timeDays - A_pc.nominalArrival.hohmannArrivalDays;

targetRadiusPatched   = A_pc.nominalArrival.periapsisRadiusKm;
targetAltitudePatched = A_pc.nominalArrival.periapsisAltitudeKm;

nbodyClosestRadius   = A_nb.closestApproach.distanceKm;
nbodyClosestAltitude = A_nb.closestApproach.altitudeKm;

radiusErrorVsNominal_km   = nbodyClosestRadius   - targetRadiusPatched;
altitudeErrorVsNominal_km = nbodyClosestAltitude - targetAltitudePatched;

% TRAJECTORY DEVIATION OVER COMMON TIME INTERVAL
tEnd = min(tNB(end), tPC(end));
maskNB = tNB <= tEnd;
maskPC = tPC <= tEnd;

tCommon = tNB(maskNB);

rPC_i = interp1(tPC(maskPC), rPC(maskPC,:), tCommon, 'linear');
vPC_i = interp1(tPC(maskPC), vPC(maskPC,:), tCommon, 'linear');

deltaR = vecnorm(rNB(maskNB,:) - rPC_i, 2, 2);
deltaV = vecnorm(vNB(maskNB,:) - vPC_i, 2, 2);

[maxDeltaR, idxMaxDeltaR] = max(deltaR);
[maxDeltaV, idxMaxDeltaV] = max(deltaV);

timeMaxDeltaR_days = tCommon(idxMaxDeltaR)/86400;
timeMaxDeltaV_days = tCommon(idxMaxDeltaV)/86400;

rmsDeltaR = sqrt(mean(deltaR.^2));
rmsDeltaV = sqrt(mean(deltaV.^2));

% Normalised trajectory deviation by heliocentric radius of n-body solution
rNBmag_common = vecnorm(rNB(maskNB,:), 2, 2);
deltaR_norm = deltaR ./ rNBmag_common;
maxDeltaR_norm = max(deltaR_norm);
rmsDeltaR_norm = sqrt(mean(deltaR_norm.^2));

% PHASE ERROR (HELIOCENTRIC ANGLE)
thetaNB = unwrap(atan2(rNB(maskNB,2), rNB(maskNB,1)));
thetaPC = unwrap(atan2(rPC_i(:,2),    rPC_i(:,1)));
deltaTheta = thetaNB - thetaPC;

[maxAbsDeltaTheta, idxMaxTheta] = max(abs(deltaTheta));
timeMaxTheta_days = tCommon(idxMaxTheta)/86400;

% METRICS AT NOMINAL HOHMANN ARRIVAL TIME
tNom = H_pc.tT;

rNB_nom = interp1(tNB, rNB, tNom, 'linear');
vNB_nom = interp1(tNB, vNB, tNom, 'linear');
rPC_nom = interp1(tPC, rPC, tNom, 'linear');
vPC_nom = interp1(tPC, vPC, tNom, 'linear');

rVenus_nom = circular_orbit_pos(R_V_pc, nV_pc, tNom, thetaV0)';
vVenus_nom = circular_orbit_vel(R_V_pc, nV_pc, tNom, thetaV0)';

deltaR_nom_km   = norm(rNB_nom - rPC_nom);
deltaV_nom_kmps = norm(vNB_nom - vPC_nom);

distVenusNB_nom_km = norm(rNB_nom - rVenus_nom);
distVenusPC_nom_km = norm(rPC_nom - rVenus_nom);

relVelVenusNB_nom_kmps = norm(vNB_nom - vVenus_nom);
relVelVenusPC_nom_kmps = norm(vPC_nom - vVenus_nom);

deltaTheta_nom_rad = atan2(rNB_nom(2), rNB_nom(1)) - atan2(rPC_nom(2), rPC_nom(1));

% TRAJECTORY DEVIATION AT ENCOUNTER-RELEVANT TIMES
tCA_PC = tPC(idxMinPC);
tCA_NB = tNB(idxMinNB);

rNB_at_PCCA = interp1(tNB, rNB, tCA_PC, 'linear');
vNB_at_PCCA = interp1(tNB, vNB, tCA_PC, 'linear');
rPC_at_PCCA = interp1(tPC, rPC, tCA_PC, 'linear');
vPC_at_PCCA = interp1(tPC, vPC, tCA_PC, 'linear');

rNB_at_NBCA = interp1(tNB, rNB, tCA_NB, 'linear');
vNB_at_NBCA = interp1(tNB, vNB, tCA_NB, 'linear');
rPC_at_NBCA = interp1(tPC, rPC, tCA_NB, 'linear');
vPC_at_NBCA = interp1(tPC, vPC, tCA_NB, 'linear');

deltaR_at_PCCA_km   = norm(rNB_at_PCCA - rPC_at_PCCA);
deltaV_at_PCCA_kmps = norm(vNB_at_PCCA - vPC_at_PCCA);

deltaR_at_NBCA_km   = norm(rNB_at_NBCA - rPC_at_NBCA);
deltaV_at_NBCA_kmps = norm(vNB_at_NBCA - vPC_at_NBCA);

% RUNTIME METRICS
runtimeRatio_NB_over_PC = A_nb.runtime / A_pc.runtime;

% PRINT SUMMARY
fprintf('====================================================\n');
fprintf('COMPARISON: EARTH-VENUS PATCHED-CONIC VS N-BODY\n');
fprintf('====================================================\n\n');

fprintf('[Closest-approach comparison]\n');
fprintf('Patched closest-approach time          = %.6f days\n', A_pc.closestApproach.timeDays);
fprintf('N-body closest-approach time           = %.6f days\n', A_nb.closestApproach.timeDays);
fprintf('Time difference (NB - PC)              = %.6f days\n\n', timeDiffClosest_days);

fprintf('Patched closest-approach distance      = %.6f km\n', A_pc.closestApproach.distanceKm);
fprintf('N-body closest-approach distance       = %.6f km\n', A_nb.closestApproach.distanceKm);
fprintf('Distance difference (NB - PC)          = %.6f km\n\n', distanceDiffClosest_km);

fprintf('Patched closest-approach altitude      = %.6f km\n', A_pc.closestApproach.altitudeKm);
fprintf('N-body closest-approach altitude       = %.6f km\n', A_nb.closestApproach.altitudeKm);
fprintf('Altitude difference (NB - PC)          = %.6f km\n\n', altitudeDiffClosest_km);

fprintf('Patched rel. speed at closest approach = %.6f km/s\n', A_pc.closestApproach.relSpeedKmps);
fprintf('N-body rel. speed at closest approach  = %.6f km/s\n', A_nb.closestApproach.relSpeedKmps);
fprintf('Rel. speed difference (NB - PC)        = %.6f km/s\n\n', relSpeedDiffClosest_kmps);

fprintf('[Nominal patched arrival target vs n-body realised outcome]\n');
fprintf('Patched nominal arrival time           = %.6f days\n', A_pc.nominalArrival.hohmannArrivalDays);
fprintf('N-body closest-approach time           = %.6f days\n', A_nb.closestApproach.timeDays);
fprintf('Timing shift vs nominal                = %.6f days\n\n', timingShiftVsNominal_days);

fprintf('Patched intended periapsis radius      = %.6f km\n', targetRadiusPatched);
fprintf('Patched intended periapsis altitude    = %.6f km\n', targetAltitudePatched);
fprintf('N-body closest-approach distance       = %.6f km\n', nbodyClosestRadius);
fprintf('N-body closest-approach altitude       = %.6f km\n', nbodyClosestAltitude);
fprintf('Radius error vs nominal patched target = %.6f km\n', radiusErrorVsNominal_km);
fprintf('Altitude error vs nominal patched tgt  = %.6f km\n\n', altitudeErrorVsNominal_km);

fprintf('[State deviation over common interval]\n');
fprintf('Max ||r_NB - r_PC||                    = %.6f km\n', maxDeltaR);
fprintf('Time of max ||r_NB - r_PC||            = %.6f days\n', timeMaxDeltaR_days);
fprintf('RMS ||r_NB - r_PC||                    = %.6f km\n', rmsDeltaR);
fprintf('Max normalised ||r_NB-r_PC||           = %.6e\n', maxDeltaR_norm);
fprintf('RMS normalised ||r_NB-r_PC||           = %.6e\n\n', rmsDeltaR_norm);

fprintf('Max ||v_NB - v_PC||                    = %.6f km/s\n', maxDeltaV);
fprintf('Time of max ||v_NB - v_PC||            = %.6f days\n', timeMaxDeltaV_days);
fprintf('RMS ||v_NB - v_PC||                    = %.6f km/s\n\n', rmsDeltaV);

fprintf('[Phase error]\n');
fprintf('Max |delta theta|                      = %.6e rad\n', maxAbsDeltaTheta);
fprintf('Time of max |delta theta|              = %.6f days\n\n', timeMaxTheta_days);

fprintf('[Metrics at nominal Hohmann arrival time]\n');
fprintf('Nominal arrival time                   = %.6f days\n', tNom/86400);
fprintf('||r_NB - r_PC|| at nominal arrival     = %.6f km\n', deltaR_nom_km);
fprintf('||v_NB - v_PC|| at nominal arrival     = %.6f km/s\n', deltaV_nom_kmps);
fprintf('NB distance to Venus at nominal time   = %.6f km\n', distVenusNB_nom_km);
fprintf('PC distance to Venus at nominal time   = %.6f km\n', distVenusPC_nom_km);
fprintf('NB rel speed to Venus at nominal time  = %.6f km/s\n', relVelVenusNB_nom_kmps);
fprintf('PC rel speed to Venus at nominal time  = %.6f km/s\n', relVelVenusPC_nom_kmps);
fprintf('Phase error at nominal arrival         = %.6e rad\n\n', deltaTheta_nom_rad);

fprintf('[State deviation at encounter-relevant epochs]\n');
fprintf('||r_NB-r_PC|| at PC closest approach   = %.6f km\n', deltaR_at_PCCA_km);
fprintf('||v_NB-v_PC|| at PC closest approach   = %.6f km/s\n', deltaV_at_PCCA_kmps);
fprintf('||r_NB-r_PC|| at NB closest approach   = %.6f km\n', deltaR_at_NBCA_km);
fprintf('||v_NB-v_PC|| at NB closest approach   = %.6f km/s\n\n', deltaV_at_NBCA_kmps);

fprintf('[Runtime]\n');
fprintf('Patched-conic runtime                  = %.6f s\n', A_pc.runtime);
fprintf('Restricted n-body runtime              = %.6f s\n', A_nb.runtime);
fprintf('Runtime ratio (NB / PC)                = %.6f\n', runtimeRatio_NB_over_PC);

% PLOTTING

% Plot 1: Heliocentric trajectory comparison
figure;
plot(rPC(:,1), rPC(:,2), 'LineWidth', 1.4); hold on;
plot(rNB(:,1), rNB(:,2), 'LineWidth', 1.4);
plot(rEarthPC(:,1), rEarthPC(:,2), '--', 'LineWidth', 1.1);
plot(rVenusPC(:,1), rVenusPC(:,2), '--', 'LineWidth', 1.1);
plot(0,0,'k.','MarkerSize',18);
plot(rPC(idxMinPC,1), rPC(idxMinPC,2), 'o', 'MarkerSize', 8, 'LineWidth', 1.4);
plot(rNB(idxMinNB,1), rNB(idxMinNB,2), 's', 'MarkerSize', 8, 'LineWidth', 1.4);
axis equal; grid on;
xlabel('x (km)');
ylabel('y (km)');
title('Earth-Venus heliocentric trajectory comparison');
legend('Patched-conic spacecraft','Restricted n-body spacecraft', ...
       'Earth orbit','Venus orbit','Sun', ...
       'PC closest approach','NB closest approach', ...
       'Location','best');

% Plot 2: Spacecraft-Venus distance vs time
figure;
plot(tPC/86400, distVenusPC, 'LineWidth', 1.4); hold on;
plot(tNB/86400, distVenusNB, 'LineWidth', 1.4);
plot(A_pc.closestApproach.timeDays, A_pc.closestApproach.distanceKm, 'o', 'MarkerSize', 8, 'LineWidth', 1.4);
plot(A_nb.closestApproach.timeDays, A_nb.closestApproach.distanceKm, 's', 'MarkerSize', 8, 'LineWidth', 1.4);
xline(tNom/86400, '--');
grid on;
xlabel('Time since SOI exit (days)');
ylabel('Spacecraft-Venus distance (km)');
title('Spacecraft-Venus distance: patched-conic vs restricted n-body');
legend('Patched-conic','Restricted n-body','PC closest approach','NB closest approach', ...
       'Nominal Hohmann arrival','Location','best');

% Plot 3: Trajectory deviation
figure;
plot(tCommon/86400, deltaR, 'LineWidth', 1.5);
grid on;
xlabel('Time (days)');
ylabel('||r_{NB} - r_{PC}|| (km)');
title('Trajectory deviation: restricted n-body vs patched-conic');

% Plot 4: Velocity deviation
figure;
plot(tCommon/86400, deltaV, 'LineWidth', 1.5);
grid on;
xlabel('Time (days)');
ylabel('||v_{NB} - v_{PC}|| (km/s)');
title('Velocity deviation: restricted n-body vs patched-conic');

% Plot 5: Phase error
figure;
plot(tCommon/86400, deltaTheta, 'LineWidth', 1.5);
grid on;
xlabel('Time (days)');
ylabel('\Delta\theta (rad)');
title('Heliocentric phase error: restricted n-body vs patched-conic');

% Plot 6: Normalised trajectory deviation
figure;
plot(tCommon/86400, deltaR_norm, 'LineWidth', 1.5);
grid on;
xlabel('Time (days)');
ylabel('||r_{NB} - r_{PC}|| / ||r_{NB}||');
title('Normalised trajectory deviation');

        % Plot 7: Zoom near closest approach
windowDays = 10;
figure;
plot(tPC/86400, distVenusPC, 'LineWidth', 1.4); hold on;
plot(tNB/86400, distVenusNB, 'LineWidth', 1.4);
xline(tNom/86400, '--');
xlim([min([A_pc.closestApproach.timeDays, A_nb.closestApproach.timeDays]) - windowDays, ...
      max([A_pc.closestApproach.timeDays, A_nb.closestApproach.timeDays]) + windowDays]);
grid on;
xlabel('Time since SOI exit (days)');
ylabel('Spacecraft-Venus distance (km)');
title('Zoomed Venus encounter comparison');
legend('Patched-conic','Restricted n-body','Nominal arrival','Location','best');

% Plot 8: Runtime comparison
figure;
bar([A_pc.runtime, A_nb.runtime]);
grid on;
set(gca, 'XTickLabel', {'Patched-conic','Restricted n-body'});
ylabel('Runtime (s)');
title('Runtime comparison: Earth-Venus case');