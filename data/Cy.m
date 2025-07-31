function Cy_t = Cy_calculation(Alpha, Beta, del_lef, del_r, del_a, b, V, p, r)
    
    function [alpha,beta,Cx] = read_data(filename)
        data = readmatrix(filename, 'Delimiter',',');
        alpha = data(:,1);
        beta = data(:,2);
        Cx = data(:,3);
    end

    function [alpha,Cx] = read_data1(filename)
        data = readmatrix(filename, 'Delimiter',',');
        alpha = data(:,1);
        Cx = data(:,2);
    end

    function Cy = Cy_alpha_beta(Alpha, Beta)
        [alpha,beta,Cy] = read_data("data\aero_coeffs\Cy\Cy(alpha,beta).txt");
        F = scatteredInterpolant(alpha, beta, Cy, 'linear');
        Cy = F(Alpha,Beta);
    end


    %fprintf('%.4f', Cy_alpha_beta(20,20));

    function Cy = Cy_lef(Alpha,Beta)
        [alpha,beta,Cy] = read_data("data\aero_coeffs\Cy\C_y_lef(alpha,beta).txt");
        F = scatteredInterpolant(alpha, beta, Cy, 'linear');
        Cy = F(Alpha, Beta);
    end
    
    %fprintf('%.4f', Cy_lef(20,20));

    function Cy = Cy_del_a_20(Alpha,Beta)
        [alpha, beta, Cy] = read_data("data\aero_coeffs\Cy\Cy_del_a=20(alpha,beta).txt");
        F = scatteredInterpolant(alpha, beta, Cy);
        Cy = F(Alpha, Beta);
    end

    %fprintf('%.4f', Cy_del_a_20(20,20));

    function Cy = del_Cy_a_20_lef(Alpha, Beta)
        [alpha, beta, Cy] = read_data("data\aero_coeffs\Cy\Cy_del_a=20,lef.txt");
        F = scatteredInterpolant(alpha, beta, Cy);
        Cy = F(Alpha, Beta);
    end
    
    %fprintf('%.4f',del_Cy_a_20_lef(20,20));

    function Cy = del_Cy_del_r_30(Alpha, Beta)
        [alpha, beta, Cy] = read_data("data\aero_coeffs\Cy\Cy,r_del_r=30(alpha,beta).txt");
        F = scatteredInterpolant(alpha, beta, Cy);
        Cy = F(Alpha, Beta);
    end

    %fprintf("%.4f", del_Cy_del_r_30(20,20));

    function Cy = Cy_r(Alpha)
        [alpha, Cy] = read_data1("data\aero_coeffs\Cy\C_y_r(alpha).txt");
        Cy = interp1(alpha, Cy, Alpha, 'linear');
    end

    %fprintf("%.4f", Cy_r(20));

    function Cy = del_Cy_r_lef(Alpha)
        [alpha, Cy] = read_data1("data\aero_coeffs\Cy\delta(C_y_r,lef).txt");
        Cy = interp1(alpha, Cy, Alpha, 'linear');
    end
    
    %fprintf("%.4f",del_Cy_r_lef(20));

    function Cy = Cy_p(Alpha)
        [alpha, Cy] = read_data1("data\aero_coeffs\Cy\C_y_p.txt");
        Cy = interp1(alpha, Cy, Alpha, 'linear');
    end

    %fprintf("%.4f",Cy_p(20));

    function Cy = del_Cy_p_lef(Alpha)
        [alpha, Cy] = read_data1("data\aero_coeffs\Cy\delta(C_y_p,lef).txt");
        Cy = interp1(alpha, Cy, Alpha, 'linear');
    end

    %fprintf("%.4f",del_Cy_p_lef(20));

    
    del_Cy_lef = Cy_lef(Alpha, Beta) - Cy_alpha_beta(Alpha,Beta);
    del_Cy_del_a_20 = Cy_del_a_20(Alpha,Beta) - Cy_alpha_beta(Alpha, Beta);


    Cy_t = Cy_alpha_beta(Alpha, Beta) + (del_Cy_lef * (1 - (del_lef/25))) + ...
            ((del_Cy_del_a_20 + del_Cy_a_20_lef(Alpha, Beta) * (1 - (del_lef/25))) * (del_a/20)) + ...
            (del_Cy_del_r_30(Alpha, Beta) * (del_r/30)) + ...
            (b/2*V)*(Cy_r(Alpha) + (del_Cy_r_lef(Alpha)*(1-(del_lef/25))))*r + ...
            (Cy_p(Alpha) + (del_Cy_p_lef(Alpha)*(1 - (del_lef/25)))) * p;

end


fprintf('%.4f', Cy_calculation(20,20, 0.4, 0.5, 0.1, 25, 150, 1000, 2000));