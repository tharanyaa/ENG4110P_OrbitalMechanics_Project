% Week 1: Two-body Hohmann transfer baseline
clear; clc; close all;

% --- Inputs (Earth example; units: km, s) ---
mu  = 3.986004418e5;        % km^3/s^2 (Earth)
Re  = 6378.1363;            % km
h1  = 300;                  % km altitude initial
h2  = 35786;                % km altitude target (GEO-ish)

r1 = Re + h1;
r2 = Re + h2;

% --- Hohmann transfer calculations ---
vc1 = sqrt(mu/r1);
vc2 = sqrt(mu/r2);

at  = (r1 + r2)/2;

vt1 = sqrt(mu*(2/r1 - 1/at));
vt2 = sqrt(mu*(2/r2 - 1/at));

dv1 = vt1 - vc1;
dv2 = vc2 - vt2;
dvTot = dv1 + dv2;

tTrans = pi*sqrt(at^3/mu);  % seconds (time of flight)

fprintf("Hohmann transfer results:\n");
fprintf("dv1   = %.4f km/s\n", dv1);
fprintf("dv2   = %.4f km/s\n", dv2);
fprintf("dvTot = %.4f km/s\n", dvTot);
fprintf("tTrans= %.2f hours\n", tTrans/3600);

% --- Plot geometry (2D) ---
theta = linspace(0,2*pi,400);
rCirc1 = r1;
rCirc2 = r2;

% Transfer ellipse in polar form with periapsis at r1
e  = (r2 - r1)/(r2 + r1);
p  = at*(1 - e^2);
thetaTrans = linspace(0, pi, 400);
rTrans = p ./ (1 + e*cos(thetaTrans));

% Convert to Cartesian
x1 = rCirc1*cos(theta); y1 = rCirc1*sin(theta);
x2 = rCirc2*cos(theta); y2 = rCirc2*sin(theta);
xT = rTrans.*cos(thetaTrans); yT = rTrans.*sin(thetaTrans);

figure;
plot(x1,y1,'LineWidth',1.5); hold on;
plot(x2,y2,'LineWidth',1.5);
plot(xT,yT,'LineWidth',2);
plot(0,0,'k.','MarkerSize',18);
axis equal; grid on;
xlabel('x (km)'); ylabel('y (km)');
legend('Initial circular orbit','Target circular orbit','Hohmann transfer','Earth');
title('Week 1: Two-body Hohmann baseline');