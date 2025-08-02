% properties og F-16

function P = F16_prop()
m = 9286; % kg
Ixx = 12875; %kg-m^2
Iyy = 75674; %kg-m^2
Izz = 85552; %kg-m^2
Ixz = 1331; %kg-m^2

b = 9.144; %m (wing span)
S = 27.87; %m^2 (wing area)
c = 3.45; %m (Mean aerodynamic chord)

x_cg_ref  = 0.35*c;
x_cg = 0.38*c;

P = [m; Ixx; Iyy; Izz; Ixz; b; S; c; x_cg_ref; x_cg];

end
