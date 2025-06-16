% Cx

function Cx_t = Cx_calculation(Alpha, Beta, Del_h, del_lef, del_sb, c, q, V)

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

    function Cx = Cx_lef(Alpha,Beta)
        [alpha, beta, Cx] = read_data('C:\Users\deepa\NDI_controller_for_a_fighter_aircraft\data\aero_coeffs\Cx\Cx,lef(alpha,beta).txt');
        F = scatteredInterpolant(alpha, beta, Cx, "linear");
        Cx = F(Alpha, Beta);
    end
    
    %fprintf('%.4f', Cx_lef(20,20));
    
    function Cx = Cx_alpha_beta_del_h_0(Alpha, Beta)
        [alpha, beta, Cx] = read_data('C:\Users\deepa\NDI_controller_for_a_fighter_aircraft\data\aero_coeffs\Cx\Cx(alpha,beta,del_h=0).txt');
        F = scatteredInterpolant(alpha, beta, Cx, "linear");
        Cx = F(Alpha, Beta);
    end

    %fprintf('%.4f', Cx_alpha_beta_del_h_0(20, 20));

    function Cx = Cx_alpha_beta_del_h(Alpha,Beta,Del_h)
        data1 = readmatrix('C:\Users\deepa\NDI_controller_for_a_fighter_aircraft\data\aero_coeffs\Cx\Cx(alpha,beta,del_h=-25).txt','Delimiter',',');
        data2 = readmatrix('C:\Users\deepa\NDI_controller_for_a_fighter_aircraft\data\aero_coeffs\Cx\Cx(alpha,beta,del_h=-10).txt','Delimiter',',');
        data3 = readmatrix('C:\Users\deepa\NDI_controller_for_a_fighter_aircraft\data\aero_coeffs\Cx\Cx(alpha,beta,del_h=0).txt','Delimiter',',');
        data4 = readmatrix('C:\Users\deepa\NDI_controller_for_a_fighter_aircraft\data\aero_coeffs\Cx\Cx(alpha,beta,del_h=10).txt','Delimiter',',');
        data5 = readmatrix('C:\Users\deepa\NDI_controller_for_a_fighter_aircraft\data\aero_coeffs\Cx\Cx(alpha,beta,del_h=25).txt','Delimiter',',');
        
        alpha = [data1(:,1); data2(:,1); data3(:,1); data4(:,1); data5(:,1)];
        beta = [data1(:,2); data2(:,2); data3(:,2); data4(:,2); data5(:,2)];
        del_h = [-25*ones(size(data1(:,1))); -10*ones(size(data2(:,1))); zeros(size(data3(:,1))); 10*ones(size(data4(:,1))); 25*ones(size(data5(:,1)))];
        C_x = [data1(:,3); data2(:,3); data3(:,3); data4(:,3); data5(:,3)];

        F = scatteredInterpolant(alpha, beta, del_h, C_x, 'linear');
        Cx = F(Alpha, Beta, Del_h);
    end

    %fprintf('%.4f', Cx_alpha_beta_del_h(20,20,25));

    function Cx = del_Cx_sb(Alpha)
        [alpha, Cx] = read_data1("C:\Users\deepa\NDI_controller_for_a_fighter_aircraft\data\aero_coeffs\Cx\delCx,sb(alpha).txt");
        Cx = interp1(alpha, Cx, Alpha, 'linear');
    end

    %fprintf('%.4f', del_Cx_sb(20));

    function Cx = Cx_q(Alpha)
        [alpha, Cx] = read_data1("C:\Users\deepa\NDI_controller_for_a_fighter_aircraft\data\aero_coeffs\Cx\Cxq(alpha).txt");
        Cx = interp1(alpha, Cx, Alpha, 'linear');
    end
    %fprintf('%.4f', Cx_q(20));

    function Cx = del_Cx_q_lef(Alpha)
        [alpha, Cx] = read_data1("C:\Users\deepa\NDI_controller_for_a_fighter_aircraft\data\aero_coeffs\Cx\delCxq,lef (alpha).txt");
        Cx = interp1(alpha, Cx, Alpha, 'linear');
    end

    %fprintf('%.4f',  del_Cx_q_lef(20));


    del_Cx_lef = Cx_lef(Alpha, Beta) - Cx_alpha_beta_del_h_0(Alpha,Beta);
    Cx_t = Cx_alpha_beta_del_h(Alpha, Beta, Del_h) + ...
        del_Cx_lef * (1 - (del_lef/25)) + ...
        del_Cx_sb(Alpha) * (del_sb/60) + ...
        (c*q/2*V) * (Cx_q(Alpha) + del_Cx_q_lef(Alpha)*(1-(del_lef/25)));

end

%fprintf('%.4f', Cx_calculation(20, 20, 25, 12, 0.4, 25, 40, 100));