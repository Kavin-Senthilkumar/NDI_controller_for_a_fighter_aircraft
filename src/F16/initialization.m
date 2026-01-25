% Initializing parameters

% mass = 9295.4128; % kg
Ixx = 12875; % kg-m^2
Iyy = 75674; % kg-m^2
Izz = 85552; % kg-m^2
Ixz = 1331; % kg-m^2

I = [Ixx 0 0; 0 Iyy 0; Ixz 0 Izz];
c_bar = 3.45; % m

b = 9.144; % m

S = 27.87; % m^2
g = 9.81; % m/s^2

xcgr = 0.35*c_bar; % m
xcg = 0.3*c_bar; % m

V_reference_mps = 300; % mps
rho_reference_kgpm3 = 1.225; % kgpm3
nv = 0.7;
nrho = 1.0;

Tmax = 105700; % N
xf_m = -0.3;
zf_m = 0.5;
alfaf_deg = -2.5; % degrees