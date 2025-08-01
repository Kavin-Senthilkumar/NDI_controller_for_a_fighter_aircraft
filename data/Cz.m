function Cz_t = Cz_calculation(Alpha, Beta, del_h, del_lef, del_sb, c, q, V)

    function [alpha, beta, Cz] = read_data(filename)
        data = readmatrix(filename,'Delimiter',',');
        alpha = data(:,1);
        beta = data(:,2);
        Cz = data(:,3);
    end
    
    function [alpha,Cz] = read_data1(filename)
        data = readmatrix(filename,'Delimiter',',');
        alpha = data(:,1);
        Cz = data(:,2);
    end

    function Cz = Cz_lef(Alpha, Beta)
        [alpha, beta, Cz] = read_data('data\aero_coeffs\Cz\Cz,lef(alpha,beta)');
        F = scatteredInterpolant(alpha, beta, Cz, "linear");
        Cz = F(Alpha, Beta);
    end

    %fprintf('%.4f', Cz_lef(45,30));

    function Cz = Cz_alpha_beta_del_h_0(Alpha, Beta)
        [alpha, beta, Cz] = read_data('data\aero_coeffs\Cz\Cz(alpha,beta,del_h=0)');
        F = scatteredInterpolant(alpha, beta, Cz, "linear");
        Cz = F(Alpha, Beta);
    end
    
    %fprintf('%.4f', Cz_alpha_beta_del_h_0(45,30));

    function Cz = Cz_alpha_beta_del_h(Alpha, Beta, Del_h)
       data1 = readmatrix('data\aero_coeffs\Cz\Cz(alpha,beta,del_h=-25)','Delimiter',',');
       data2 = readmatrix('data\aero_coeffs\Cz\Cz(alpha,beta,del_h=-10)','Delimiter',',');
       data3 = readmatrix('data\aero_coeffs\Cz\Cz(alpha,beta,del_h=0)','Delimiter',',');
       data4 = readmatrix('data\aero_coeffs\Cz\Cz(alpha,beta,del_h=10)','Delimiter',',');
       data5 = readmatrix('data\aero_coeffs\Cz\Cz(alpha,beta,del_h=25)','Delimiter',',');

       alpha = [data1(:,1); data2(:,1); data3(:,1); data4(:,1); data5(:,1)];
       beta = [data1(:,2); data2(:,2); data3(:,2); data4(:,2); data5(:,2)];
       del_h = [-25*ones(size(data1(:,1))); -10*ones(size(data2(:,1))); zeros(size(data3(:,1))); 10*ones(size(data4(:,1))); 25*ones(size(data5(:,1)))];
       C_z = [data1(:,3); data2(:,3); data3(:,3); data4(:,3); data5(:,3)];

       F = scatteredInterpolant(alpha, beta, del_h, C_z, 'linear');
       Cz = F(Alpha, Beta, Del_h);
    end

    %fprintf('%.4f', Cz_alpha_beta_del_h(45, 20, 25));

    function Cz = del_Cz_sb_alpha(Alpha)
        [alpha, Cz] = read_data1('data\aero_coeffs\Cz\delCz,sb(alpha)');
        Cz = interp1(alpha, Cz, Alpha, "linear");
    end
%fprintf('%.4f', del_Cz_sb_alpha(25));

    function Cz = Czq_alpha(Alpha)
        [alpha, Cz] = read_data1('data\aero_coeffs\Cz\Czq(alpha).txt');
        Cz = interp1(alpha, Cz, Alpha, 'linear');
    end
%fprintf('%.4f', Czq_alpha(25));

    function Cz = del_Czq_lef(Alpha)
        [alpha, Cz] = read_data1('data\aero_coeffs\Cz\delCzq,lef(alpha)');
        Cz = interp1(alpha, Cz, Alpha, 'linear');
    end

%fprintf('%.4f', del_Czq_lef(25));

del_Cz_lef = Cz_lef(Alpha, Beta) - Cz_alpha_beta_del_h_0(Alpha, Beta);


Cz_t = Cz_alpha_beta_del_h(Alpha, Beta, del_h) + ...
    del_Cz_lef*(1 - (del_lef/25)) + ...
    del_Cz_sb_alpha(Alpha) * (del_sb/60) + ...
    ((c*q)/(2*V))*(Czq_alpha(Alpha) + (del_Czq_lef(Alpha)*(1-(del_lef/25))));
   

end

