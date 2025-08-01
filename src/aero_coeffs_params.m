
function aero_params = aero_coeffs_params()
    Alpha = 5.0;
    Beta = 2.0;
    del_a = 1.0;
    del_r = 0.5;
    del_h = -2.0;
    del_sb = 0.0;
    del_lef = 10.0;
    p = 0.02;
    q = 0.01;
    r = 0.015;
    c = 2.5;
    b = 15.0;
    V = 100;
    x_cg = 0.3;
    x_cg_ref = 0.25;

    aero_params = [Alpha; Beta; del_a; del_r; del_h; del_sb; del_lef ; p; q; r; c; b; V; x_cg; x_cg_ref];
end