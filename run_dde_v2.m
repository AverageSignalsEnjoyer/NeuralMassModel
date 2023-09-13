%%

clear all
close all
rng('default')
rng(1)

%% Specify network parameters 

% ---- Excitatory and inhibitory synaptic gain ----
He = 3.25; % mV
Hi = 29.3; % mV

% ---- P.Time constants ----
Te = 10/1000; % s
Ti = 15/1000; % s

% ---- Average number of synaptic contacts ----
g1 = 50;
g2 = 40;
g3 = 12;
g4 = g3;

% ---- Sigmoid func parameters ----
e0 = 2.5;  % Hz   - 1/2 of max firing rate
v0 = 0;    % mV   - potential when 50% firing rate achieved
r  = 0.56; % 1/mV - slope of sigmoid function at v0

% ---- Input parameters ----
C     = 1000;
sigma = 0.05;

% ---- Simulation parameters ----
fs   = 100;
dur  = 3;
tVec = 1/fs:1/fs:dur;

% ---- Network parameters ----
dim = [5 5];   % spatial dimension of the network
N = prod(dim); % total number of columns
ns = 8;        % number of state parameters per column

% ---- Ephaptic coupling parameters ----
d  = 2/1000;   % m     - length of pyramidal cell(?);
gn = 10e-14;   % mA/mV - input-to-current (Mercadal 2022)
R  = 0.4;      % S/m   - conductivity (Sanchez-Todo 2022) 
ds = 5*d;      % m     - unit of distance in the grid

%% Build connectivity

% ---- Define the inputs into the network ----
u = zeros(N, length(tVec));
u(1, 10:100:end) = 1;
% u(N, 10:100:end) = 1;
% u = randn([N, length(tVec)])*sigma; % Random noise input

% ---- Coupling and delay matrices ----
Cf = zeros(N, N); % forward
Cb = zeros(N, N); % backward
Cl = zeros(N, N); % lateral
D  = zeros(N, N); % ms - delay between pairs of columns 
E  = zeros(N, N); % ephaptic coupling factor (Nunez 5.6)

% % Linear chain
% for i = 1:N
%     for j = 1:N
%         if i - j == 1
%             Cf(i, j) = 10;
%             D(i, j) = 0.01;
%         end
%     end
% end

% 2D Grid - Neighbors - NW -> SE
for i = 1:N
    for j = 1:N
        
        [m, n] = ind2sub(dim, i);
        [p, q] = ind2sub(dim, j);

        if p-m == 1 && n == q
            Cf(j, i) = 20;
            Cb(i, j) = 5; % uncomment for feedback
            D(i, j) = 0.01;
            D(j, i) = 0.01;
        end

        if q-n == 1 && p == m
            Cf(j, i) = 20;
            Cb(i, j) = 5; % uncomment for feedback
            D(i, j) = 0.01;
            D(j, i) = 0.01;
        end

        if q-n == 1 && p-m == 1
            Cf(j, i) = 20;
            Cb(i, j) = 5; % uncomment for feedback
            D(i, j) = 0.01;
            D(j, i) = 0.01;
        end

        tempDist = sqrt((ds*(p-m))^2 + (ds*(q-n))^2);
        rad = sqrt((d/2)^2 + tempDist^2);

        if i ~= j
            E(j, i) = d^2/(8*pi*R*rad^3);
        end

    end
end

% ---- Conduction delay matrices ----

% State equation matrices
A  = zeros(N*ns, N*ns);     % A*x      - Normal state parameters
As = zeros(N*ns, N*ns);     % A*S(x)   - Sigmoid transformed parameters
Ay = zeros(N*ns, N);        % A*S(y)   - Within-column pyramidal cell
Au = zeros(N*ns, N);        % A*u      - Extrinsic inputs
Ad = zeros(N, N, N*ns, N);  % A*S(y-d) - Inter-column pyramidal cell
Ae = zeros(N*ns, N);        % A*x      - Ephaptic coupling

for n = 1:N

    i = (n-1)*ns;

    A(i+1, i+4) = 1;
    A(i+2, i+5) = 1;
    A(i+3, i+6) = 1;
    A(i+7, i+8) = 1;

    A(i+4, i+4) = -2/Te;
    A(i+4, i+1) = -1/Te^2;
    Au(i+4, n)  = He/Te*C;
    Ay(i+4, n)  = He/Te*g1;

    A(i+5, i+5)  = -2/Te;
    A(i+5, i+2)  = -1/Te^2;
    As(i+5, i+1) = He/Te*g2;

    A(i+6, i+6)  = -2/Ti;
    A(i+6, i+3)  = -1/Ti^2;
    As(i+6, i+7) = Hi/Ti*g4;

    A(i+8, i+8) = -2/Te;
    A(i+8, i+7) = -1/Te^2;
    Ay(i+8, n)  = He/Te*g3;

    % ---- Inter-region interactions ----
    for m = 1:N
        if m ~= n
            Ad(n, m, i+4, m) = He/Te*(Cf(n, m) + Cl(n, m));
            Ad(n, m, i+5, m) = He/Te*(Cb(n, m) + Cl(n, m));
            Ad(n, m, i+8, m) = He/Te*(Cb(n, m) + Cl(n, m));
        end
    end

    % ---- Ephaptic coupling ----
    % for m = 1:N
    %     if m ~= n
    %         idx = (m-1)*ns;
    %         Ae(i+5, m) =  He/Te*E(n, m);
    %         Ae(i+6, m) = -Hi/Ti*E(n, m);
    %     end
    % end

end

%% Run diff eq

tSpan = [1/fs tVec(end)];
init  = zeros(1, N*ns);
opts  = odeset('MaxStep', 0.01);

lags = unique(D);
lags(lags <= 0) = [];

% ---- Packpage params to pass to dde ----
P.A  = A;   
P.As = As;
P.Ay = Ay;
P.Au = Au;
P.Ad   = Ad;
P.Ae   = Ae;
P.D    = D;
P.lags = lags;
P.tVec = tVec;
P.e0   = e0;
P.v0   = v0;
P.r    = r;
P.N    = N;
P.ns   = ns;
P.u    = u;

sol = dde23(@(t,x,z) fn_dde_v2(t, x, z, P), lags, init, tSpan, opts);

y_  = resample(sol.y', sol.x, fs)';
y = zeros(N, size(y_, 2));

% Calculate pyramidal cell depolarization
for i = 1:N
    idx = (i-1)*ns;
    y(i, :) = y_(idx+2, :) - y_(idx+3, :);
end

figure
hold on
plot(tVec, y')

figure
nT = 25;
temp = y(:, 1:floor(size(y, 2)/nT):size(y, 2));
for i = 1:nT
    subplot(5, 5, i)
    imagesc(reshape(temp(:, i), dim))
end

save('y_fb', 'y', 'sol', 'dim', 'fs', 'tVec')
