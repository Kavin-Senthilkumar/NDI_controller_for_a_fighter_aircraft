% Simulating F-16 aircraft

% Physical parameters
global CGref CGact;
m = 8936;       % Mass (kg)
Ixx = 12896;    % Moment of inertia about x-axis (kg.m^2)
Iyy = 75800;    % Moment of inertia about y-axis (kg.m^2)
Izz = 85672;    % Moment of inertia about z-axis (kg.m^2)
Ixz = 1332;     % Product of inertia (kg.m^2)
S = 27.87;      % Wing area (m^2)
b = 9.96;       % Wing span (m)
c = 3.45;       % Mean aerodynamic chord (m)
CGref = 1.2075; % Reference center of gravity (m)
CGact = 1.311;  % Actual center of gravity (m)
Hx = 216.9;     % Angular momentum (kg.m^2/s) - not used in this model
g = 9.81;       % Gravitational acceleration (m/s^2)

% ADC Function: Computes air data (dynamic pressure and Mach number)
function air_data = ADC(V_t, alt)
    rho0 = 1.225;          % Sea-level air density (kg/m^3)
    H = 10400;             % Scale height (m)
    rho = rho0 * exp(-alt / H); % Air density at altitude
    qbar = 0.5 * rho * V_t^2;   % Dynamic pressure (Pa)
    speed_of_sound = 340;  % Speed of sound (m/s)
    mach = V_t / speed_of_sound; % Mach number
    air_data.qbar = qbar;
    air_data.mach = mach;
end

% ADCE Function: Computes aerodynamic coefficients
function aero_coeff = ADCE(qbar, V_t, c, b, p, q, r, alpha, beta, ELV, RDR, AIL)
    global CGref CGact;
    alpha_deg = rad2deg(alpha); % Convert angle of attack to degrees
    beta_deg = rad2deg(beta);   % Convert sideslip angle to degrees
    
    % Baseline aerodynamic coefficients
    Cx0 = -0.032 - (0.08 * alpha_deg) - (0.02 * ELV); % Longitudinal force coefficient
    Cy0 = -(0.57 * beta_deg) + (0.05 * RDR / 30);     % Lateral force coefficient
    Cz0 = -0.34 - (1.8 * alpha_deg) - (0.2 * ELV);    % Normal force coefficient
    Cl0 = -(0.02 * beta_deg) + (0.15 * AIL / 20) + (0.02 * RDR / 30); % Roll moment coefficient
    Cm0 = 0.05 - (0.21 * alpha_deg) - (0.5 * ELV);    % Pitch moment coefficient
    Cn0 = (0.057 * beta_deg) + (0.01 * AIL / 20) + (0.19 * RDR / 30); % Yaw moment coefficient
    
    % Total aerodynamic coefficients with dynamic effects
    Cx_t = Cx0; % No dynamic terms assumed for Cx
    Cz_t = Cz0; % No dynamic terms assumed for Cz
    Cl_t = Cl0 + ((b / (2 * V_t)) * (0.08 * r - 0.41 * p)); % Roll moment with rate terms
    Cm_t = Cm0 + ((c / (2 * V_t)) * q * -1.2) + (Cz0 * (CGref - CGact)); % Pitch moment with CG shift
    Cn_t = Cn0 + ((b / (2 * V_t)) * (-0.15 * r - 0.01 * p)) - (Cy0 * (CGref - CGact) * c / b); % Yaw moment with CG shift
    
    aero_coeff = [Cx_t Cz_t Cl_t Cm_t Cn_t]; % Return coefficients as a vector
end

% dof Function: Defines the 12-state equations of motion
function dy = dof(t, y, c, b, AIL, RDR, ELV, m, S, Ixx, Iyy, Izz, Ixz)
    global CGref CGact;
    % State variables
    V_t = y(1);    % Total velocity (m/s)
    alpha = deg2rad(y(2)); % Angle of attack (rad)
    beta = deg2rad(y(3));  % Sideslip angle (rad)
    p = y(4);      % Roll rate (rad/s)
    q = y(5);      % Pitch rate (rad/s)
    r = y(6);      % Yaw rate (rad/s)
    phi = deg2rad(y(7));   % Roll angle (rad)
    theta = deg2rad(y(8)); % Pitch angle (rad)
    psi = deg2rad(y(9));   % Yaw angle (rad)
    xpos = y(10);  % X position (m)
    ypos = y(11);  % Y position (m)
    zpos = y(12);  % Z position (m, positive downward)
    g = 9.81;
    % Compute air data
    AD = ADC(V_t, -zpos); % Altitude is negative zpos
    qbar = AD.qbar;       % Dynamic pressure
    
    % Compute aerodynamic coefficients
    aero_coeff = ADCE(qbar, V_t, c, b, p, q, r, alpha, beta, ELV, RDR, AIL);
    Cx = aero_coeff(1); % Longitudinal force coefficient
    Cz = aero_coeff(2); % Normal force coefficient
    Cl = aero_coeff(3); % Roll moment coefficient
    Cm = aero_coeff(4); % Pitch moment coefficient
    Cn = aero_coeff(5); % Yaw moment coefficient
    
    % Velocity components in body frame
    u = V_t * cos(alpha) * cos(beta); % X-axis velocity (m/s)
    v = V_t * sin(beta);              % Y-axis velocity (m/s)
    w = V_t * sin(alpha) * cos(beta); % Z-axis velocity (m/s)
    
    % Thrust (simplified as zero for this model)
    T = 0;         % Thrust force (N)
    sigmaT = 0;    % Thrust angle (rad)
    
    % Force equations (acceleration in body frame)
    udot = -(q * w) + (r * v) - (g * sin(theta)) + (qbar * S * Cx) / m + (T * cos(sigmaT)) / m;
    vdot = -(r * u) + (p * w) + (g * cos(theta) * sin(phi));
    wdot = -(p * v) + (q * u) + (g * cos(theta) * cos(phi)) + (qbar * S * Cz) / m - (T * sin(sigmaT)) / m;
    
    % Moment equations (angular acceleration)
    pdot = (qbar * S * b * (Izz * Cl + Ixz * Cn) - q * r * (Ixz^2 + Izz^2 + Iyy * Izz) + p * q * Ixz * (Ixx - Iyy + Izz)) / ((Ixx * Izz) - (Ixz^2));
    qdot = (qbar * S * c * Cm - (p^2 - r^2) * Ixz + p * r * (Izz - Ixx)) / Iyy;
    rdot = (qbar * S * b * (Ixx * Cn + Ixz * Cl) - q * r * (Ixx - Iyy + Izz) + p * q * (Ixz^2 + Ixx^2 - Ixx * Iyy)) / (Ixx * Izz - Ixz^2);
    
    % Euler angle rates
    phidot = p + q * tan(theta) * sin(phi) + r * tan(theta) * cos(phi);
    thetadot = q * cos(phi) - r * sin(phi);
    psidot = (r * cos(phi) + q * sin(phi)) * sec(theta);
    
    % Position rates in inertial frame
    xdot = u * cos(theta) * cos(psi) + v * (sin(phi) * sin(theta) * cos(psi) - cos(phi) * sin(psi)) + w * (cos(phi) * sin(theta) * cos(psi) + sin(phi) * sin(psi));
    ydot = u * cos(theta) * sin(psi) + v * (sin(phi) * sin(theta) * sin(psi) + cos(phi) * cos(psi)) + w * (cos(phi) * sin(theta) * sin(psi) - sin(phi) * cos(psi));
    zdot = -u * sin(theta) + v * sin(phi) * cos(theta) + w * cos(phi) * cos(theta);
    
    % State derivative vector
    dy = [udot; vdot; wdot; pdot; qdot; rdot; phidot; thetadot; psidot; xdot; ydot; zdot];
end

% Main Script: Solves the equations of motion and plots results
% Initial conditions: [V_t, alpha_deg, beta_deg, p, q, r, phi_deg, theta_deg, psi_deg, x, y, z]
y0 = [100; 4; 1; 0.1; 0.05; 0.08; 5; 0; 10; 500; 0; -200];
tspan = [0 10]; % Time span (seconds)
ELV = 0.0;       % Elevator deflection (degrees)
RDR = 0.0;       % Rudder deflection (degrees)
AIL = 0.0;       % Aileron deflection (degrees)

% Solve the differential equations using ode45
[t, y] = ode23(@(t, y) dof(t, y, c, b, AIL, RDR, ELV, m, S, Ixx, Iyy, Izz, Ixz), tspan, y0);

% Plotting results
% 3D Position Plot
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