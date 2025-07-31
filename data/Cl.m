function Cl_t = Cl_calculation(Alpha, Beta, del_lef, del_a, del_r, del_h, r, p, b, V)

% Read data for 3D interpolation (alpha, beta)
function [alpha, beta, Cl] = read_data(filename)
    data = readmatrix(filename, 'Delimiter', ',');
    alpha = data(:,1);
    beta = data(:,2);
    Cl = data(:,3);
end

% Read data for 1D interpolation (alpha)
function [alpha, Cl] = read_data_1(filename)
    data = readmatrix(filename, 'Delimiter', ',');
    alpha = data(:,1);
    Cl = data(:,2);
end

% Cl_lef_alpha_beta
function Cl = Cl_lef_alpha_beta(Alpha, Beta)
    [alpha, beta, Cl] = read_data('data\aero_coeffs\Cl\C_l,lef(alpha,beta).txt');
    F = scatteredInterpolant(alpha, beta, Cl, 'linear');
    Cl = F(Alpha, Beta);
end

% Cl_alpha_beta_del_h_0
function Cl = Cl_alpha_beta_del_h_0(Alpha, Beta)
    [alpha, beta, Cl] = read_data('data\aero_coeffs\Cl\Cl(alpha,beta,del_h=0).txt');
    F = scatteredInterpolant(alpha, beta, Cl, 'linear');
    Cl = F(Alpha, Beta);
end

% Cl_del_a_20
function Cl = Cl_del_a_20(Alpha, Beta)
    [alpha, beta, Cl] = read_data('data\aero_coeffs\Cl\C_l,del_a=20(alpha,beta).txt');
    F = scatteredInterpolant(alpha, beta, Cl, 'linear');
    Cl = F(Alpha, Beta);
end

% Cl_del_a_20_lef
function Cl = Cl_del_a_20_lef(Alpha, Beta)
    [alpha, beta, Cl] = read_data('data\aero_coeffs\Cl\C_l,del_a=20,lef(alpha,beta).txt');
    F = scatteredInterpolant(alpha, beta, Cl, 'linear');
    Cl = F(Alpha, Beta);
end

% Cl_del_r_30
function Cl = Cl_del_r_30(Alpha, Beta)
    [alpha, beta, Cl] = read_data('data\aero_coeffs\Cl\C_l,del_r=30(alpha,beta).txt');
    F = scatteredInterpolant(alpha, beta, Cl, 'linear');
    Cl = F(Alpha, Beta);
end

% Cl_p
function Cl = Cl_p(Alpha)
    [alpha, Cl] = read_data_1('data\aero_coeffs\Cl\Cl_p(alpha).txt');
    Cl = interp1(alpha, Cl, Alpha, 'linear');
end

% Cl_r
function Cl = Cl_r(Alpha)
    [alpha, Cl] = read_data_1('data\aero_coeffs\Cl\C_l_r(alpha).txt');
    Cl = interp1(alpha, Cl, Alpha, 'linear');
end

% del_Cl_p_lef
function Cl = del_Cl_p_lef(Alpha)
    [alpha, Cl] = read_data_1('data\aero_coeffs\Cl\del_Cl_p,lef(alpha).txt');
    Cl = interp1(alpha, Cl, Alpha, 'linear');
end

% del_Cl_beta
function Cl = del_Cl_beta(Alpha)
    [alpha, Cl] = read_data_1('data\aero_coeffs\Cl\del_Cl_beta(alpha).txt');
    Cl = interp1(alpha, Cl, Alpha, 'linear');
end

% del_Cl_r_lef
function Cl = del_Cl_r_lef(Alpha)
    [alpha, Cl] = read_data_1('data\aero_coeffs\Cl\del_Cl_r,lef(alpha).txt');
    Cl = interp1(alpha, Cl, Alpha, 'linear');
end

% Cl_alpha_beta_del_h
function Cl = Cl_alpha_beta_del_h(Alpha, Beta, Del_h)
    data1 = readmatrix('data\aero_coeffs\Cl\Cl(alpha,beta,del_h=-25).txt', 'Delimiter', ',');
    data2 = readmatrix('data\aero_coeffs\Cl\Cl(alpha,beta,del_h=0).txt', 'Delimiter', ',');
    data3 = readmatrix('data\aero_coeffs\Cl\Cl(alpha,beta,del_h=25).txt', 'Delimiter', ',');
    
    alpha = [data1(:,1); data2(:,1); data3(:,1)];
    beta = [data1(:,2); data2(:,2); data3(:,2)];
    del_h = [-25*ones(size(data1(:,1))); zeros(size(data2(:,1))); 25*ones(size(data3(:,1)))];
    Cl = [data1(:,3); data2(:,3); data3(:,3)];
    
    F = scatteredInterpolant(alpha, beta, del_h, Cl, 'linear');
    Cl = F(Alpha, Beta, Del_h);
end

% Main calculations
del_Cl_lef = Cl_lef_alpha_beta(Alpha, Beta) - Cl_alpha_beta_del_h_0(Alpha, Beta);
del_Cl_del_a_20 = Cl_del_a_20(Alpha, Beta) - Cl_alpha_beta_del_h_0(Alpha, Beta);
del_Cl_del_a_20_lef = Cl_del_a_20_lef(Alpha, Beta) - Cl_lef_alpha_beta(Alpha, Beta) - (Cl_del_a_20(Alpha, Beta) - Cl_alpha_beta_del_h_0(Alpha, Beta));
del_Cl_del_r_30 = Cl_del_r_30(Alpha, Beta) - Cl_alpha_beta_del_h_0(Alpha, Beta);

Cl_t = Cl_alpha_beta_del_h(Alpha, Beta, del_h) + ...
       ((del_Cl_lef * (1 - (del_lef/25))) + ...
       (del_Cl_del_a_20 + del_Cl_del_a_20_lef * (1 - (del_lef/25))) * (del_a/20)) + ...
       (del_Cl_del_r_30 * (del_r/30)) + ...
       ((b/(2*V)) * (Cl_r(Alpha) + (del_Cl_r_lef(Alpha) * (1 - (del_lef/25)))) * r + ...
       (Cl_p(Alpha) + (del_Cl_p_lef(Alpha) * (1 - (del_lef/25))) * p)) + ...
       del_Cl_beta(Alpha) * Beta;

%fprintf('The total rolling moment coefficient is %f\n', Cl_t);

end


fprintf('.%4f', Cl_calculation(10,2,0.5,0.3,0.5,20, 25, 5, 6, 150));
