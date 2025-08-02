 %Control surface limits

function u = control_saturations(u1,u2,u3,u4,u5,u6)
    u1_max = deg2rad(21.5); %Aileron deflection
    u1_min = deg2rad(-21.5);

    u2_max = deg2rad(30); % Rudder deflection
    u2_min = deg2rad(-30);

    u3_max = deg2rad(25); % Elevator deflection
    u3_min = deg2rad(-25);

    u4_max = deg2rad(25); % Leading Edge flap
    u4_min = deg2rad(-2);

    u5_max = deg2rad(60); % Speed Brakes
    u5_min = deg2rad(0);

    u6_min = deg2rad(0); % Throttle
    u6_max = deg2rad(90);

    if u1>u1_max
        u1 = u1_max;
    elseif u1<u1_min
        u1 = u1_min;
    end

    if u2>u2_max
        u2 = u2_max;
    elseif u2<u2_min
        u2 = u2_min;
    end

    if u3>u3_max
        u3 = u3_max;
    elseif u3<u3_min
        u3 = u3_min;
    end
    
    if u4>u4_max
        u4 = u4_max;
    elseif u4<u4_min
        u4 = u4_min;
    end

    if u5>u5_max
        u5 = u5_max;
    elseif u5<u1_min
        u5 = u5_min;
    end

    if u6>u6_max
        u6 = u6_max;
    elseif u6<u6_min
        u6 = u6_min;
    end

    u = [u1; u2; u3; u4; u5; u6];
end