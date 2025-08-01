% Aerodynamic Model of the F-16 Aircraft
% addpath('data\')

function aero_coeffs = aero(Alpha, Beta, del_h, del_a, c, q, qbar, V, r, p, b, del_a, del_r, del_sb)


    Cl_t = Cl_calculation(Alpha, Beta, del_h, del_a, del_r, del_lef, r, p, b, V);
    Cx_t = Cx_calculation(Alpha, Beta, del_h, del_lef, del_sb, c, q, V);
    Cm_t = Cm_calculation(Alpha,Beta,del_h, del_lef, del_sb, x_cg_ref, x_cg, c, q, qbar, V);
    Cy_t = Cy_calculation(Alpha, Beta, del_lef, del_r, del_a, b, V, p, r);
    Cn_t = Cn_calculation(Alpha, Beta, del_h, del_lef, del_a, c, b, x_cg_ref, x_cg, r, V, beta);
    Cz_t = Cz_calculation(Alpha, Beta, del_h, del_lef, del_sb, c, q, V);




    Cl_t = Cl_calculation(Alpha, Beta, del_lef, del_a, del_r, del_h, r, p, b, V);
    Cx_t = Cx_calculation(Alpha, Beta, del_h, del_lef, del_sb, c, q, V);
    Cm_t = Cm_calculation(Alpha,Beta,del_h, del_lef, del_sb, x_cg_ref, x_cg, c, q, qbar, V);
    Cy_t = Cy_calculation(Alpha, Beta, del_lef, del_r, del_a, b, V, p, r);
    Cn_t = Cn_calculation(Alpha, Beta, del_h, del_lef, del_a, c, b, x_cg_ref, x_cg, r, V, beta);
    Cz_t = Cz_calculation(Alpha, Beta, del_h, del_lef, del_sb, c, q, V);
    



end


