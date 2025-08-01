



function Cm_t = Cm_calculation(Alpha,Beta,Del_h, del_lef, del_sb, x_cg_ref, x_cg, c, q, qbar, V)

    function [alpha,beta,Cm] = read_data(filename)
        data = readmatrix(filename,'Delimiter',',');
        alpha = data(:,1);
        beta = data(:,2);
        Cm = data(:,3);
    end

    function [alpha,Cm] = read_data_1(filename)
        data = readmatrix(filename, 'Delimiter', ',');
        alpha = data(:, 1);
        Cm = data(:, 2);
    end
    function Cm = C_m_lef_alpha_beta(Alpha,Beta)
        [alpha,beta,Cm] = read_data('data\aero_coeffs\Cm\C_m,lef(alpha,beta).txt');
        F = scatteredInterpolant(alpha,beta,Cm,"linear");
        Cm = F(Alpha,Beta); 
    end

%fprintf('%.4f',C_m_lef_alpha_beta(-5,20));

    function Cm = Cm_alpha_beta_del_h_0(Alpha,Beta)
        [alpha,beta,Cm] = read_data('data\aero_coeffs\Cm\Cm(alpha,beta,del_h = 0).txt');
        F = scatteredInterpolant(alpha, beta, Cm, "linear");
        Cm = F(Alpha, Beta);
    end
 %fprintf('%.4f',Cm_alpha_beta_del_h_0(-5,20));

    function Cm = del_C_m_sb_alpha(Alpha)
        % Implementation for del_C_m_sb_alpha function
        [alpha, Cm] = read_data_1('data\aero_coeffs\Cm\del_Cm_sb(alpha).txt');
        Cm = interp1(alpha, Cm, Alpha, 'linear', 'extrap');
    end
 %fprintf('%.4f',del_C_m_sb_alpha(20));
 
    function Cm = Cmq_alpha(Alpha)
        [alpha, Cm] = read_data_1('data\aero_coeffs\Cm\Cm_q(alpha).txt');
        Cm = interp1(alpha, Cm, Alpha, 'linear');
    end
    
    %fprintf('%.4f',Cmq_alpha(20));

    function Cm = del_Cmq_lef_alpha(Alpha)
        [alpha, Cm] = read_data_1('data\aero_coeffs\Cm\del_Cm_q_lef(alpha).txt');
        Cm = interp1(alpha, Cm, Alpha, 'linear');
    end
    %fprintf('%.4f',del_Cmq_lef_alpha(20));
    function Cm = del_Cm_alpha(Alpha)
        [alpha, Cm] = read_data_1('data\aero_coeffs\Cm\del_Cm(alpha).txt');
        Cm = interp1(alpha, Cm, Alpha, 'linear');
    end
    %fprintf('%.4f',del_Cm_alpha(20));
    
    function Cm = del_Cm_ds_alpha_del_h(Alpha, Del_h)
        data = readmatrix('data\aero_coeffs\Cm\del_C_m,ds(alpha,del_h).txt','Delimiter',',');
        alpha = data(:,1);
        del_h = data(:,2);
        Cm = data(:,3);
        F = scatteredInterpolant(alpha, del_h, Cm, "linear");
        Cm = F(Alpha, Del_h);
    end
    %fprintf('%.4f',del_Cm_ds_alpha_del_h(20,20));

    function eta_del_h = etaDel_h(Del_h)
        data = readmatrix('data\aero_coeffs\Cm\eta_del_h(del_h).txt','Delimiter',',');
        del_h = data(:,1);
        Eta_del_h = data(:,2);
        eta_del_h = interp1(del_h, Eta_del_h, Del_h, 'linear');

    end

    %fprintf('%.4f',etaDel_h(25));

    function Cm = Cm_alpha_beta_del_h(Alpha, Beta, Del_h)
        data1 = readmatrix('data\aero_coeffs\Cm\Cm(alpha,beta,del_h = -25).txt','Delimiter',',');
        data2 = readmatrix('data\aero_coeffs\Cm\Cm(alpha,beta,del_h = -10).txt','Delimiter',',');
        data3 = readmatrix('data\aero_coeffs\Cm\Cm(alpha,beta,del_h = 0).txt','Delimiter',',');
        data4 = readmatrix('data\aero_coeffs\Cm\Cm(alpha,beta,del_h=10).txt','Delimiter',',');
        data5 = readmatrix('data\aero_coeffs\Cm\Cm(alpha,beta,del_h = 25).txt','Delimiter',',');

        alpha = [data1(:,1); data2(:,1); data3(:,1); data4(:,1); data5(:,1)];
        beta = [data1(:,2); data2(:,2); data3(:,2); data4(:,2); data5(:,2)];
        del_h = [-25*ones(size(data1(:,1))); -10*ones(size(data2(:,1))); zeros(size(data3(:,1))); 10*ones(size(data4(:,1))); 25*ones(size(data5(:,1)))];
        C_m = [data1(:,3); data2(:,3); data3(:,3); data4(:,3); data5(:,3)];

        F = scatteredInterpolant(alpha, beta, del_h, C_m, 'linear');
        Cm = F(Alpha,Beta,Del_h);
    end

    %fprintf('%.4f\n',Cm_alpha_beta_del_h(20,20,25));

    del_Cm_lef = C_m_lef_alpha_beta(Alpha, Beta) - Cm_alpha_beta_del_h_0(Alpha, Beta);

    Cz_t =  Cz(Alpha, Beta, Del_h, del_lef, del_sb, c, q, V);

    Cm_t = Cm_alpha_beta_del_h(Alpha,Beta,Del_h) * etaDel_h(Del_h) + ...
        Cz_t * (x_cg_ref - x_cg) + ... %Need to import Cz_t from another file.
        del_Cm_lef(1 - (del_lef/25)) + ...
        del_C_m_sb_alpha(Alpha) * (del_sb/60) +...
        ((c*qbar)/(2*V)) * ((Cmq_alpha(Alpha) + del_Cmq_lef_alpha(Alpha))*(1-(del_lef/25))) + ...
        del_Cm_alpha(Alpha) + ...
        del_Cm_ds_alpha_del_h(Alpha, Del_h);

end
