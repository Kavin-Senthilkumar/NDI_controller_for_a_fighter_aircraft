%Simulating 6DOF equations of an aircraft without any controller

function dydt = dof(t,y)
u = y(1);v = y(2);w = y(3);p = y(4) ;q= y(5);r = y(6); phi = deg2rad(y(7)); theta = deg2rad(y(8)); psi = deg2rad(y(9)); xpos = y(10);ypos = y(11); zpos = y(12);

Cx = 0.05;
Cy = 0.02;
Cz = 0.1;

m = 1000;
g = 9.81;

sigmaT = deg2rad(10);
qbar = 250;
S= 2;
T = 5000;

Ixx = 10000;
Iyy = 20000;
Izz = 30000;
Ixz = 500;
b = 15;
c = 3;
Cl = 0.005;
Cm = -0.02;
Cn = 0.01;
%change in linear velocity - Force equations
udot = -(q*w) + (r*v) -(g*sin(theta)) +(qbar*S*Cx)/m + (T*cos(sigmaT))/m;
vdot = -(r*u) + (p*w) + (g*cos(theta)*sin(phi));
wdot = -(p*v) + (q*u) + (g*cos(theta)*cos(phi)) + (qbar*S*Cz)/m - (T*sin(sigmaT))/m;

%change in angular velocity - Moment equations
pdot = (qbar*S*b*(Izz*Cl + Ixz*Cn) - q*r*(Ixz^2 + Izz^2 + Iyy*Izz) + p*q*Ixz*(Ixx - Iyy + Izz))/((Ixx*Izz) - (Ixz^2));
qdot = (qbar*S*c*Cm - (p^2-q^2)*Ixz + p*r*(Izz - Ixx))/Iyy;
rdot = (qbar*S*b*(Ixx*Cn + Ixz*Cl) - q*r*(Ixx - Iyy + Izz) + p*q*(Ixz^2 + Ixx^2 - Ixx*Iyy))/(Ixx*Izz - Ixz^2);

%Euler angles - attitude of the aircraft
phidot = p + q*tan(theta)*sin(phi) + r*tan(theta)*cos(phi);
thetadot = (q*cos(theta))-(r*sin(phi));
psidot = (r*cos(phi)*sec(theta))+(q*sin(phi)*sec(theta));

%position vectors
xdot = u*(cos(theta)*cos(psi)) + v*(sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi)) + w*(cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi));
ydot = u*(cos(theta)*sin(psi)) + v*(sin(phi)*sin(theta)*sin(psi)) + w*(cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(phi));
zdot = -u*sin(theta) + v*sin(phi)*cos(theta) + w*cos(phi)*cos(theta);


dydt = transpose([udot vdot wdot pdot qdot rdot phidot thetadot psidot xdot ydot zdot]);
end



y0 = [50 5 -2 0.1 0.05 0.08 5 0 10 0 0 0];

tspan = [0 1000];
[t,y] = ode45(@dof,tspan,y0);


u = y(:,1);
v = y(:,2);
w = y(:,3);

p = y(:,4);
q = y(:,5);
r = y(:,6);

phi = y(:,7);
theta = y(:,8);
psi = y(:,9);

xpos = y(:,10);
ypos = y(:,11);
zpos = y(:,12);

figure;
plot3(xpos,ypos,zpos);

figure;
plot(t,u, t,v, t,w);
xlabel('Time');
ylabel('Linear Velocity');
legend('u','v','w');

figure;
plot(t,p, t,q, t,r);
xlabel('Time');
ylabel('Angular Velocity');
legend('p','q','r');

figure;
plot(t,phi, t,theta, t,psi);
xlabel('Time');
ylabel('Euler Angles');
legend('phi','theta','psi');