function [System_params, GFMInv_params] = setting_params(G_VQ)
    % network parameters
    System_params.Lg = 0.1885/(100*pi);
    System_params.Rg = 0.01; % default 0.01
    
    System_params.theta_inf = 0;
    System_params.V_inf = 1;
    
    System_params.omega0_nominal = 100*pi;

    System_params.Xg = System_params.Lg * System_params.omega0_nominal;
    
    % device parameters
    GFMInv_params.Lf = 0.012/(100*pi);
    GFMInv_params.Rf = 0.001;
    GFMInv_params.Cf = 0.044/(100*pi);
    
    GFMInv_params.p_ref = 1.0;
    GFMInv_params.v_ref = 1;
    GFMInv_params.q_ref = 0.0;
    GFMInv_params.M = 5;
    GFMInv_params.D = 10;
    GFMInv_params.Dv = 1;
    GFMInv_params.kip = 0.3;
    GFMInv_params.kii = 7.854;
    GFMInv_params.kvp = 0.1;
    GFMInv_params.kvi = 0.1;
    GFMInv_params.tau = 0.002;

    % voltage control method
    if nargin < 1
        G_VQ = GFMInv_params.Dv;
    end

    % Q-V droop transfer function: y(s) / u(s) = num(s) / den(s)
    if isnumeric(G_VQ)
        G_VQ = tf(G_VQ);
    end
    [num, den] = tfdata(G_VQ, 'v');
    [A_d, B_d, C_d, D_d] = tf2ss(num, den);
    GFMInv_params.droop_A = A_d;
    GFMInv_params.droop_B = B_d;
    GFMInv_params.droop_C = C_d;
    GFMInv_params.droop_D = D_d;
end