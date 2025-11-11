% Six DOF Differential equation

function dy = F16model(Y,U)


    % State variables
    u = Y(1); % y1
    v = Y(2); % y2
    w = Y(3); % y3
    
    p = Y(4); % y4
    q = Y(5); % y5
    r = Y(6); % y6

    phi = deg2rad(Y(7));   % Roll angle (rad)
    theta = deg2rad(Y(8)); % Pitch angle (rad)
    psi = deg2rad(Y(9)); % Yaw angle (rad)

    xpos = Y(10);  % X position (m)
    ypos = Y(11);  % Y position (m)
    zpos = Y(12);  % Z position (m, positive downward)

    del_A = U(1);
    del_R = U(2);
    del_E = U(3);
    del_LEF = U(4);
    del_SB = U(5);
    thrtl = U(6);

    V_t = sqrt(u^2 + v^2 + w^2);
    alpha = atan2(w,u); % Angle of attack (rad)
    beta = asin(v/V_t);  % Sideslip angle (rad)
    
    F16_prop = F16_properties();

    m = F16_prop(1);
    Ixx = F16_prop(2);
    Iyy = F16_prop(3);
    Izz = F16_prop(4);
    Ixz = F16_prop(5);
    b = F16_prop(6);
    S = F16_prop(7);
    c = F16_prop(8);
    x_cg_ref = F16_prop(9);
    x_cg = F16_prop(10);

    g = 9.81;
    H_e = 216.9; % Engine Angular Momentum kg-m^2/sec
    % Compute air data
    AD = ADC(V_t, -zpos); % Altitude is negative zpos
    qbar = AD.qbar;       % Dynamic pressure

    U_controlled = control_saturation(del_A, del_R, del_E, del_LEF, del_SB, thrtl);

    del_A = deg2rad(U_controlled(1));
    del_R = deg2rad(U_controlled(2));
    del_E = deg2rad(U_controlled(3));
    del_LEF = deg2rad(U_controlled(4));
    del_SB = deg2rad(U_controlled(5));
    thrtl = deg2rad(U_controlled(6));

    
    
    % Compute aerodynamic coefficients
    aero_coeff = F16_aero(alpha, beta, del_A, del_R, del_E, del_LEF, del_SB, p, q, r, c, b, V_t, x_cg, x_cg_ref);
    Cx = aero_coeff(1); % Longitudinal force coefficient
    Cy = aero_coeff(2); % Normal force coefficient
    Cz = 0.456; %aero_coeff(3);
    Cl = 0.342; %aero_coeff(4); % Roll moment coefficient
    Cm = 0.345; %aero_coeff(5); % Pitch moment coefficient
    Cn = 0.234; %aero_coeff(6); % Yaw moment coefficient
    
    % Velocity components in body frame
    
    % Thrust (simplified as zero for this model)
    T = 0;         % Thrust force (N)
    %sigmaT = 0;    % Thrust angle (rad)
    
    % Force equations (acceleration in body frame)
    udot = -(q * w) + (r * v) - (g * sin(theta)) + (qbar * S * Cx) / m + (T/m);
    vdot = -(r * u) + (p * w) + (g * cos(theta) * sin(phi)) + ((qbar*S)/m)*Cy;
    wdot = -(p * v) + (q * u) + (g * cos(theta) * cos(phi)) + (qbar * S * Cz) / m;
    
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