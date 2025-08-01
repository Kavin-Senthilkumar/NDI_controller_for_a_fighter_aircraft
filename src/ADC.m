
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
