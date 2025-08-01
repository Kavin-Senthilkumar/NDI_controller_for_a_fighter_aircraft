% Aerodynamic Model of the F-16 Aircraft
% addpath('data\')

function aero_coeffs = aero(Alpha, Beta, del_a, del_r, del_h, del_lef, del_sb, p,q,r, c,b, V, x_cg, x_cg_ref)


    
    Cx_t = Cx(Alpha, Beta, del_h, del_lef, del_sb, c, q, V);
    Cy_t = Cy(Alpha, Beta, del_lef, del_r, del_a, b, V, p, r);
    Cz_t = Cz(Alpha, Beta, del_h, del_lef, del_sb, c, q, V);
    Cl_t = Cl(Alpha, Beta, del_h, del_lef, del_a, del_r, r, p, b, V);
    Cm_t = Cm(Alpha,Beta,del_h, del_lef, del_sb, x_cg_ref, x_cg, c, q, V);
    Cn_t = Cn(Alpha, Beta, del_h, del_lef, del_r, del_a, c, b, x_cg_ref, x_cg, p, r, V);
    

    aero_coeffs = [Cx_t; Cy_t; Cz_t; Cl_t; Cm_t; Cn_t];
end

%fprintf('%.4f, %.4f, %.4f, %.4f, %.4f, %.4f',aero(5.0, 2.0, 1.0, 0.5, -2.0, 10.0, 0.0, 0.02, 0.01, 0.015, 2.5, 15.0, 100, 0.3, 0.25));
