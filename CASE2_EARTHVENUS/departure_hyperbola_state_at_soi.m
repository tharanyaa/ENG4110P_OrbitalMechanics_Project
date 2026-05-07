function [r_rel, v_rel, hyp] = departure_hyperbola_state_at_soi(muP, rp, vInf_vec, rSOI)
% departure_hyperbola_state_at_soi

% Returns the planet-centred spacecraft state at the point where the departure hyperbola reaches r = rSOI.

% Inputs:
%   muP      : planet gravitational parameter [km^3/s^2]
%   rp       : hyperbola periapsis radius [km]
%   vInf_vec : outbound hyperbolic excess velocity vector [3x1] km/s
%   rSOI     : sphere-of-influence radius [km]

% Outputs:
%   r_rel    : planet-centred position at SOI exit [3x1] km
%   v_rel    : planet-centred velocity at SOI exit [3x1] km/s
%   hyp      : struct of useful hyperbola parameters

vInf = norm(vInf_vec);
sHat = vInf_vec / vInf;              

% Hyperbola parameters
e = 1 + rp*vInf^2/muP;
h = rp * sqrt(vInf^2 + 2*muP/rp);

beta      = acos(1/e);
theta_inf = acos(-1/e);              % = pi - beta

% Solve for the true anomaly at which r = rSOI
% r = h^2/mu / (1 + e cos(theta))
p = h^2/muP;
cos_theta_soi = (p/rSOI - 1)/e;


cos_theta_soi = max(-1, min(1, cos_theta_soi));

theta_soi = acos(cos_theta_soi);     % outbound, positive angle

% Build basis in the ecliptic plane:
% outbound asymptote lies at +theta_inf from periapsis direction pHat
% so p_Hat is s_Hat rotated by -theta_inf
pHat = rotz2d(sHat, -theta_inf);
qHat = rotz2d(pHat, pi/2);

% Position and velocity in perifocal form, mapped into inertial frame
r_rel = rSOI * (cos(theta_soi)*pHat + sin(theta_soi)*qHat);
v_rel = (muP/h) * (-sin(theta_soi)*pHat + (e + cos(theta_soi))*qHat);

% Store extra outputs
hyp.e         = e;
hyp.h         = h;
hyp.beta      = beta;
hyp.theta_inf = theta_inf;
hyp.theta_soi = theta_soi;
hyp.vInf      = vInf;
end