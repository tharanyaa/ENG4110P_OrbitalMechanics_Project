function plot_orbit_and_metrics_week2(out, titleStr)
% plot_orbit_and_metrics_week2
% Week 2 diagnostics in ONE figure (5 panels):
% 1) XY trajectory
% 2) Radius vs time (absolute)
% 3) Radius deviation
% 4) Energy vs time (absolute)
% 5) Energy deviation

t = out.t;
Y = out.Y;

r = Y(:,1:3);
v = Y(:,4:6);
mu = out.params.mu;

rnorm = sqrt(sum(r.^2, 2));
eps = specific_energy_two_body(r, v, mu);


dr   = rnorm - mean(rnorm);
deps = eps - eps(1);

% ---------- One figure, multiple panels ----------
fig = figure('Name', titleStr, 'Color', 'w');

tl = tiledlayout(fig, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, titleStr, 'Interpreter', 'none');

% 1) Trajectory (XY)
nexttile(tl, 1);
plot(r(:,1), r(:,2), 'LineWidth', 1.5); hold on;
plot(0, 0, 'k.', 'MarkerSize', 14);
axis equal; grid on;
xlabel('x (km)'); ylabel('y (km)');
title('Trajectory (XY)');

% 2) Radius vs time (absolute)
nexttile(tl, 2);
plot(t/3600, rnorm, 'LineWidth', 1.5);
grid on;
xlabel('Time (hours)'); ylabel('|r| (km)');
title('Orbital radius vs time');

% 3) Radius deviation
nexttile(tl, 3);
plot(t/3600, dr, 'LineWidth', 1.5);
grid on;
xlabel('Time (hours)'); ylabel('\Delta r (km)');
title('Radius deviation');

% 4) Energy vs time (absolute)
nexttile(tl, 4);
plot(t/3600, eps, 'LineWidth', 1.5);
grid on;
xlabel('Time (hours)'); ylabel('\epsilon (km^2/s^2)');
title('Specific energy vs time');

% 5) Energy deviation (span bottom row)
nexttile(tl, 5, [1 2]);
plot(t/3600, deps, 'LineWidth', 1.5);
grid on;
xlabel('Time (hours)'); ylabel('\Delta\epsilon (km^2/s^2)');
title('Energy deviation');
end
