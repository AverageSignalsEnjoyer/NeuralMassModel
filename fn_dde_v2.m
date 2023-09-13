function xdot = fn_dde_v2(t, x, z, P)

N  = P.N;
ns = P.ns;
nl = length(P.lags);

% ---- Interpolate inputs ----
u_ = zeros(N, 1);
for n = 1:N
    u_(n) = interp1(P.tVec, P.u(n, :), t, 'nearest'); % Also interp inputs
end
u = u_;

% ---- Calculate pyramidal cell depolarization ----
y  = zeros(N, 1);
yd = zeros(N, nl);
for n = 1:N
    idx = (n-1)*ns;
    y(n) = x(idx+2) - x(idx+3);
    for d = 1:nl
        yd(n, d) = z(idx+2, d) - z(idx+3, d);
    end
end

% ---- Calculate diff ----
xdot = P.A*x + (P.Ae+P.Ay)*S(y, P.e0, P.r) + P.Au*u + P.As*S(x, P.e0, P.r);

% ---- Calculate coupling effects ----
for n = 1:length(u)
    for m = 1:length(u)
        tempLag = P.D(n, m);
        if tempLag ~= 0
            iLag = P.lags == tempLag;
            xdot = xdot + squeeze(P.Ad(n, m, :, :)) * S(yd(:, iLag), P.e0, P.r);
        end
    end
end

end

function out = S(v, e0, r)

out = 2*e0./(1+exp(-r*v)) - e0;

end