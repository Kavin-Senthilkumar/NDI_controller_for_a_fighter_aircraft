y0 = [100; 4; 1; 0.1; 0.05; 0.08; 5; 0; 10; 500; 0; -200];

Alpha = 5.0;
Beta = 2.0;
AIL = 1.0;
RDR = 0.5;
ELV = -2.0;
SB  = 0.0;
LEF = 10.0;
m = 8936;       % Mass (kg)
Ixx = 12896;    % Moment of inertia about x-axis (kg.m^2)
Iyy = 75800;    % Moment of inertia about y-axis (kg.m^2)
Izz = 85672;    % Moment of inertia about z-axis (kg.m^2)
Ixz = 1332;     % Product of inertia (kg.m^2)
S = 27.87;      % Wing area (m^2)
b = 9.96;       % Wing span (m)
c = 3.45;       % Mean aerodynamic chord (m)
x_cg_ref = 1.2075; % Reference center of gravity (m)
x_cg = 1.311;  % Actual center of gravity (m)
Hx = 216.9;     % Angular momentum (kg.m^2/s) - not used in this model
g = 9.81; 

tspan = [0 10];

[t, y] = ode23(@(t,y) F16_deq(t, y, c, b, AIL, RDR, ELV,LEF, SB, m, S, Ixx, Iyy, Izz, Ixz, x_cg, x_cg_ref), tspan, y0);


figure;
plot3(y(:,10), y(:,11), y(:,12));
xlabel('X Position (m)');
ylabel('Y Position (m)');
zlabel('Z Position (m)');
title('F-16 3D Trajectory');
grid on;

% Linear Velocities Plot
figure;
plot(t, y(:,1), t, y(:,2), t, y(:,3));
xlabel('Time (s)');
ylabel('Linear Velocity (m/s)');
title('Linear Velocities');
legend('u', 'v', 'w');
grid on;

% Angular Velocities Plot
figure;
plot(t, y(:,4), t, y(:,5), t, y(:,6));
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');
title('Angular Velocities');
legend('p', 'q', 'r');
grid on;

% Euler Angles Plot
figure;
plot(t, y(:,7), t, y(:,8), t, y(:,9));
xlabel('Time (s)');
ylabel('Euler Angles (deg)');
title('Euler Angles');
legend('phi', 'theta', 'psi');
grid on;