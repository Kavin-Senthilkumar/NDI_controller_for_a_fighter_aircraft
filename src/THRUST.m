function THRUST = THRUST(POW, ALT, MACH)

Mach = [0.2, 0.4, 0.6, 0.8, 1.0];

altitudes = [0, 3048, 6096, 9144, 12192, 15240];

T_idle = [
    2824, 1890, 3069, 4492, 5916, 7562;
    267, 111, 1535, 3358, 5026, 6783;
    -4537, -3158, -1334, 1557, 4048, 6049;
    -12010, -8451, -5782, -1099, 2669, 4893;
    -16013, -6227, -2647, -1521, -890, 3114;
    ];

T_mil = [
    56401, 40699, 28080, 17970, 10987, 6227;
    56089, 41420, 29401, 19082, 11565, 6939;
    56223, 43764, 31536, 20728, 12632, 7384;
    55111, 45263, 34472, 23663, 14456, 8585;
    51953, 43804, 35806, 27133, 16902, 10275;
    ];

T_max = [
    95276, 69834, 49929, 32573, 19727, 11565;
    100970, 74993, 54488, 36269, 22240, 12610;
    107820, 84112, 61204, 41300, 25354, 14300;
    115959, 93742, 71057, 49440, 30513, 17570;
    128485, 103723, 81398, 59977, 38440, 22494;
    ];



% Perform 2D interpolation using interp2 without normalization
T_idle_interp = interp2(altitudes, Mach, T_idle, ALT, MACH, 'linear');
T_mil_interp = interp2(altitudes, Mach, T_mil, ALT, MACH, 'linear');
T_max_interp = interp2(altitudes, Mach, T_max, ALT, MACH, 'linear');

% Calculate thrust based on power
if POW < 50
    THRUST = T_idle_interp + (T_mil_interp - T_idle_interp) * POW * 0.02;
else
    THRUST = T_mil_interp + (T_max_interp - T_mil_interp) * (POW - 50) * 0.02;
end


end