function x_equi = calc_equilibrium(GFMInv_params, System_params, opt)
    % Virtual inductor parameter (default: not employed)
    if nargin <= 2
        opt.VIenable = false;
        opt.VIxv = 0.1;
    end

    vod_init = GFMInv_params.v_ref;
    voq_init = 0.0;
    xvd_init = 0.0;
    xvq_init = 0.0;
    iLd_init = GFMInv_params.p_ref / 1.5 / vod_init; % Approximate from P = 1.5 * vd * id
    iLq_init = GFMInv_params.q_ref / 1.5 / vod_init; % Approximate from Q = -1.5 * vd * iq (if voq=0)
    xid_init = 0.0;
    xiq_init = 0.0;
    theta_init = 0.0;
    omega_init = 0.0;
    i_LineD = 0.0;
    i_LineQ = 0.0;
    
    n_droop = size(GFMInv_params.droop_A, 1);
    x_droop_init = zeros(n_droop, 1);
        
    x_init = [vod_init; voq_init; xvd_init; xvq_init; iLd_init; iLq_init; xid_init; xiq_init; theta_init; omega_init; i_LineD; i_LineQ; x_droop_init];

    x_equi = fsolve(@(x) calc_dynamics(0, x, GFMInv_params, System_params, opt), x_init);
end