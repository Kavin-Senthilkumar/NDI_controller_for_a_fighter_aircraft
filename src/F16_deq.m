% Six DOF Differential equation

function dy = dof(~, y, c, b, AIL, RDR, ELV, m, S, Ixx, Iyy, Izz, Ixz)
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
    H_e = 216.9; % Engine Angular Momentum kg-m^2/sec
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
    pdot = (qbar * S * b * (Izz * Cl + Ixz * Cn) - q * r * (Ixz^2 + Izz^2 + Iyy * Izz) + p * q * Ixz * (Ixx - Iyy + Izz)) / ((Ixx * Izz) - (Ixz^2)) - (H_e * r);
    qdot = (qbar * S * c * Cm - (p^2 - r^2) * Ixz + p * r * (Izz - Ixx)) / Iyy;
    rdot = ((qbar * S * b * (Ixx * Cn + Ixz * Cl) - q * r * (Ixx - Iyy + Izz) + p * q * (Ixz^2 + Ixx^2 - Ixx * Iyy)) / (Ixx * Izz - Ixz^2)) + (H_e * q);
    
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