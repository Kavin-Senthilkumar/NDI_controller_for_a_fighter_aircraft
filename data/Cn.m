% Cn

function Cn_t = Cn_calculation(Alpha, Beta, Del_h, del_lef, del_a, c, b, x_cg_ref, x_cg, r, V, beta)

    function [alpha,beta,Cn] = read_data(filename)
        data = readmatrix(filename, 'Delimiter',',');
        alpha = data(:,1);
        beta = data(:,2);
        Cn = data(:,3);
    end

    function [alpha,Cn] = read_data1(filename)
        data = readmatrix(filename, 'Delimiter',',');
        alpha = data(:,1);
        Cn = data(:,2);
    end

    function Cn = Cn_lef(Alpha,Beta)
        [alpha,beta,Cn] = read_data('data\aero_coeffs\Cn\C_n,lef(alpha,beta).txt');
        F = scatteredInterpolant(alpha,beta,Cn,"linear");
        Cn = F(Alpha,Beta);
    end
    %fprintf('%.4f',Cn_lef(20,20));

    function Cn = Cn_alpha_beta_del_h_0(Alpha,Beta)
        [alpha,beta,Cn] = read_data('data\aero_coeffs\Cn\Cn(alpha,beta,del_h = 0).txt');
        F = scatteredInterpolant(alpha,beta,Cn,"linear");
        Cn = F(Alpha,Beta);
    end
    %fprintf('%.4f',Cn_alpha_beta_del_h_0(20,20));

    function Cn = Cn_del_a_20(Alpha,Beta)
        [alpha,beta,Cn] = read_data('data\aero_coeffs\Cn\Cn,del_a=20(alpha,beta).txt');
        F = scatteredInterpolant(alpha,beta,Cn,"linear");
        Cn = F(Alpha,Beta);
    end

    %fprintf('%.4f',Cn_del_a_20(20,20));
    
    function Cn = Cn_del_a_20_lef(Alpha,Beta)
        [alpha,beta,Cn] = read_data('data\aero_coeffs\Cn\C_n,del_a=20,lef(alpha,beta).txt');
        F = scatteredInterpolant(alpha,beta,Cn,"linear");
        Cn = F(Alpha,Beta);
    end
    
    %fprintf('%.4f',Cn_del_a_20_lef(20,20));

    function Cn = Cn_del_r_30(Alpha,Beta)
        [alpha,beta,Cn] = read_data('data\aero_coeffs\Cn\Cn,del_r=30(alpha,beta).txt');
        F = scatteredInterpolant(alpha,beta,Cn,"linear");
        Cn = F(Alpha,Beta);
    end

    %fprintf('%.4f',Cn_del_r_30(20,20));

    function Cn = Cn_alpha_beta_del_h(Alpha,Beta,Del_h)
        data1 = readmatrix('data\aero_coeffs\Cn\Cn(alpha,beta,del_h=-25).txt', 'Delimiter', ',');
        data2 = readmatrix('data\aero_coeffs\Cn\Cn(alpha,beta,del_h = 0).txt', 'Delimiter', ',');
        data3 = readmatrix('data\aero_coeffs\Cn\Cn(alpha,beta,del_h = 25).txt', 'Delimiter', ',');

        alpha = [data1(:,1);data2(:,1);data3(:,1)];
        beta = [data1(:,2);data2(:,2);data3(:,2)];
        del_h = [-25*ones(size(data1(:,1)));zeros(size(data2(:,1)));25*ones(size(data3(:,1)))];
        Cn = [data1(:,3);data2(:,3);data3(:,3)];
        F = scatteredInterpolant(alpha,beta,del_h,Cn);
        Cn = F(Alpha, Beta, Del_h);
    end
    
    %fprintf('%.4f',Cn_alpha_beta_del_h(20,20,25));

    function Cn = Cn_r(Alpha)
        [alpha,Cn] = read_data1('data\aero_coeffs\Cn\Cn_r(alpha).txt');
        Cn = interp1(alpha,Cn,Alpha,'linear');
    end
    %fprintf('%.4f',Cn_r(30));

    function Cn = del_Cn_r_lef(Alpha)
        [alpha,Cn] = read_data1('data\aero_coeffs\Cn\del_Cn_r,lef(alpha).txt');
        Cn = interp1(alpha,Cn,Alpha,'linear');
    end

    %fprintf('%.4f',del_Cn_r_lef(30));

    function Cn = Cn_p(Alpha)
        [alpha, Cn] = read_data1('data\aero_coeffs\Cn\Cn_p(alpha).txt');
        Cn = interp1(alpha, Cn, Alpha, 'linear');
    end

    %fprintf('%.4f',Cn_p(30));

    function Cn = del_Cn_p_lef(Alpha)
        [alpha, Cn] = read_data1('data\aero_coeffs\Cn\del_Cn_p,lef(alpha).txt');
        Cn = interp1(alpha, Cn, Alpha, 'linear');
    end

   %fprintf('%.4f',del_Cn_p_lef(30));

    function Cn = del_Cn_beta(Alpha)
        [alpha, Cn] = read_data1('data\aero_coeffs\Cn\del_Cn_beta(alpha).txt');
        Cn = interp1(alpha,Cn, Alpha, "linear");
    end

    
    del_Cn_lef = Cn_lef(Alpha,Beta) - Cn_alpha_beta_del_h_0(Alpha,Beta);
    del_Cn_del_a_20 = Cn_del_a_20(Alpha,Beta) - Cn_alpha_beta_del_h_0(Alpha,Beta);
    del_Cn_del_a_20_lef = Cn_del_a_20_lef(Alpha,Beta) - Cn_lef(Alpha,Beta) - (Cn_del_a_20(Alpha,Beta) - Cn_alpha_beta_del_h_0(Alpha,Beta));
    del_Cn_del_r_30 = Cn_del_r_30(Alpha,Beta) - Cn_alpha_beta_del_h_0(Alpha,Beta);



    Cn_t = Cn_alpha_beta_del_h(Alpha,Beta,Del_h) + ...
        del_Cn_lef*(1-(del_lef/25)) + ...
        Cy_t * (x_cg_ref - x_cg) * (c/b) + ...
        ((del_Cn_del_a_20 + ((del_Cn_del_a_20_lef)*(1-(del_lef/25)))) * (del_a/20)) + ...
        (del_Cn_del_r_30 * (del_r/30)) + ...
        (b/2*V)*(Cn_r(Alpha) + (del_Cn_r_lef(Alpha)*(1-(del_lef/25)))*r) + ...
        (Cn_p(Alpha) + (del_Cn_p_lef(Alpha) * (1-(del_lef/25)))*p) + ...
        del_Cn_beta(Alpha) * Beta ;

end
