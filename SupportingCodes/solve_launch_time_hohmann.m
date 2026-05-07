
function tL = solve_launch_time_hohmann(thetaInt0, thetaTar0, n1, n2, tT, deltaAim)
% solve_launch_time_hohmann with angular offset for targeting

if nargin < 6   %no of input arguments
    deltaAim = -0.2;
end
% Modified for aiming offset
rhs0 = (thetaTar0 + deltaAim - thetaInt0) + n2*tT - pi;

den = (n1 - n2);
if abs(den) < 1e-12
    error("n1 and n2 are too close; this formula needs different handling.");
end

k0 = ceil((-rhs0)/(2*pi));

candidates = [];
for k = (k0-3):(k0+3)
    t = (rhs0 + 2*pi*k)/den;
    if t > 0
        candidates(end+1) = t; %
    end
end

if isempty(candidates)
    for k = k0:(k0+1000)  %If nearby search failed, keep trying larger integer numbers until a positive launch time is found
        t = (rhs0 + 2*pi*k)/den;
        if t > 0
            tL = t;
            return;
        end
    end
    error("No positive launch time found within search range.");
end

tL = min(candidates);
end