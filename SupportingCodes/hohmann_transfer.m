function H = hohmann_transfer(r1, r2, mu)
% hohmann_transfer  Computes classic coplanar Hohmann transfer quantities.
% Inputs:
%   r1, r2 : radii from central body (km)
%   mu     : gravitational parameter (km^3/s^2)
% Output struct H:
%   H.aT      transfer semi-major axis (km)
%   H.tT      transfer time (s) (half period of transfer ellipse)
%   H.vc1     circular speed at r1 (km/s)
%   H.vc2     circular speed at r2 (km/s)
%   H.vpT     speed at perigee of transfer ellipse (at r1) (km/s)
%   H.vaT     speed at apogee of transfer ellipse (at r2) (km/s)
%   H.dv1     departure burn magnitude (km/s)
%   H.dv2     arrival burn magnitude (km/s)

H.aT  = 0.5*(r1 + r2);
H.tT  = pi*sqrt(H.aT^3/mu);

H.vc1 = sqrt(mu/r1);
H.vc2 = sqrt(mu/r2);

% vis-viva on transfer ellipse
H.vpT = sqrt(mu*(2/r1 - 1/H.aT));
H.vaT = sqrt(mu*(2/r2 - 1/H.aT));

H.dv1 = H.vpT - H.vc1;
H.dv2 = H.vc2 - H.vaT;
end
