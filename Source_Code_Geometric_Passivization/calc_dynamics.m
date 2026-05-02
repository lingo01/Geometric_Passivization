function dx = calc_dynamics(t, x, GFMInv_params, System_params, opt)
    %% preparations
    if nargin <= 4
        opt.disturb = 'p_ref';
        opt.VIenable = false;
        opt.VIxv = 0.1;
    end

    % droop parameters
    A_droop = GFMInv_params.droop_A;
    B_droop = GFMInv_params.droop_B;
    C_droop = GFMInv_params.droop_C;
    D_droop = GFMInv_params.droop_D;
    n_droop = size(A_droop, 1);

    % define the order of state variants
    dx = zeros(12 + n_droop, 1); % GFM + line dynamics have 12 state variables + droop states
   
    % bias of GFM and system states
    v_od = x(1);
    v_oq = x(2);
    x_vd = x(3); % voltage loop d-axis integrator state
    x_vq = x(4); % voltage loop q-axis integrator state
    i_Ld = x(5); % filter inductor current d-axis
    i_Lq = x(6); % filter inductor current q-axis
    x_id = x(7); % current loop d-axis integrator state
    x_iq = x(8); % current loop q-axis integrator state
    theta_gfm = x(9); % GFM internal angle (angle of d-axis relative to global D-axis)
    omega_gfm = x(10);% GFM rotor angular frequency deviation (rad/s, p.u. of omega_sync if M is in p.u. time)
    i_LineD = x(11); % line current D-axis component
    i_LineQ = x(12); % line current Q-axis component
    if n_droop > 0
        x_droop = x(13:12+n_droop); % droop states
    else
        x_droop = [];
    end

    %% paramter deinitions
    % GFM dynamic parameters
    p_ref = GFMInv_params.p_ref;
    v_ref = GFMInv_params.v_ref;
    q_ref = GFMInv_params.q_ref;
    Cf    = GFMInv_params.Cf;
    Lf    = GFMInv_params.Lf;
    kip   = GFMInv_params.kip;
    kii   = GFMInv_params.kii;
    kvp   = GFMInv_params.kvp;
    kvi   = GFMInv_params.kvi;
    dv    = GFMInv_params.Dv; % Q-V droop coefficient
    d_damping = GFMInv_params.D; % P-omega damping coefficient
    Rf = GFMInv_params.Rf; % filter resistance (if included in model)
    omega0_nominal = System_params.omega0_nominal; % synchronous angular frequency (e.g., 2*pi*50 or 1 p.u.)
    M     = GFMInv_params.M; % inertia

    % system dynamic paramters
    Lg = System_params.Lg;
    rg = System_params.Rg; % line resistance

    theta_inf = System_params.theta_inf;
    V_inf = System_params.V_inf;

    % set disturbance
    if t >= 1
        if strcmp(opt.disturb, 'p_ref')
            p_ref = p_ref * 1.01;
        elseif strcmp(opt.disturb, 'q_ref')
            q_ref = q_ref * 1.01;
        elseif strcmp(opt.disturb, 'v_ref')
            v_ref = v_ref * 1.01;
        else
            
        end
    end

    %% framework transformation
    % GFM output voltage in (D,Q) frameworks
    v_gfm_D = v_od * cos(theta_gfm) - v_oq * sin(theta_gfm);
    v_gfm_Q = v_od * sin(theta_gfm) + v_oq * cos(theta_gfm);

    % GFM output current in (d,q) frameworks
    i_od = i_LineD * cos(theta_gfm) + i_LineQ * sin(theta_gfm);
    i_oq = -i_LineD * sin(theta_gfm) + i_LineQ * cos(theta_gfm);

    % GFM output voltage in (D,Q) frameworks
    V_inf_D = V_inf * cos(theta_inf);
    V_inf_Q = V_inf * sin(theta_inf);

    % active and reactive power comuptation
    p_gfm = 1.5 * (v_od * i_od + v_oq * i_oq);
    q_gfm = 1.5 * (v_oq * i_od - v_od * i_oq); % Note: q = vq*id - vd*iq
    
    %% GFM internal dynamics and line dynamics
    % Swing Equation
    dx(10) = (p_ref - p_gfm - d_damping * (omega_gfm)) / M; % d(omega_gfm)/dt
    dx(9)  = omega0_nominal * omega_gfm; % d(theta_gfm)/dt, theta_gfm is absolute angle

    % Q-V droop states calculation
    u_droop = q_gfm - q_ref;
    if n_droop > 0
        dx(13:12+n_droop) = A_droop * x_droop + B_droop * u_droop;
        y_droop = C_droop * x_droop + D_droop * u_droop;
    else
        y_droop = D_droop * u_droop;
    end

    % Voltage Control Loop
    if ~opt.VIenable % if not employing the virtual inductor control
        vod_ref = v_ref - y_droop; % Q-V droop
        voq_ref = 0; % usually q-axis voltage reference is 0
    else
        vod_ref = v_ref - y_droop + opt.VIxv * i_Lq; % Q-V droop
        voq_ref = 0 - opt.VIxv * i_Ld; % usually q-axis voltage reference is 0
    end

    dx(3) = vod_ref - v_od; % d(x_vd)/dt
    dx(4) = voq_ref - v_oq; % d(x_vq)/dt

    % Current reference calculation (from voltage loop output)
    iLd_ref = kvp * (vod_ref - v_od) + kvi * x_vd - omega0_nominal * Cf * v_oq + i_od;
    iLq_ref = kvp * (voq_ref - v_oq) + kvi * x_vq + omega0_nominal * Cf * v_od + i_oq;

    iLd_ref = calc_limiter(iLd_ref, [-1.5, 1.5]);
    iLq_ref = calc_limiter(iLq_ref, [-1.5, 1.5]);

    % Current control loop
    dx(7) = iLd_ref - i_Ld; % d(x_id)/dt
    dx(8) = iLq_ref - i_Lq; % d(x_iq)/dt

    % PWM voltage calculation (current loop output)
    v_id_pwm = kip * (iLd_ref - i_Ld) + kii * x_id + v_od - omega0_nominal * Lf * i_Lq;
    v_iq_pwm = kip * (iLq_ref - i_Lq) + kii * x_iq + v_oq + omega0_nominal * Lf * i_Ld;

    v_id_pwm = calc_limiter(v_id_pwm, [-1.2, 1.2]);
    v_iq_pwm = calc_limiter(v_iq_pwm, [-1.2, 1.2]);

    % LC filter dynamics (inductor current)
    dx(5) = (v_id_pwm - v_od - Rf * i_Ld + omega0_nominal * (1+omega_gfm) * Lf * i_Lq) / Lf;
    dx(6) = (v_iq_pwm - v_oq - Rf * i_Lq - omega0_nominal * (1+omega_gfm) * Lf * i_Ld) / Lf;

    % LC filter dynamics (capacitor voltage)
    dx(1) = (i_Ld - i_od + omega0_nominal * (1+omega_gfm) * Cf * v_oq) / Cf;
    dx(2) = (i_Lq - i_oq - omega0_nominal * (1+omega_gfm) * Cf * v_od) / Cf;

    % Line dynamics (global DQ frame)
    dx(11) = (v_gfm_D - V_inf_D + omega0_nominal*Lg*i_LineQ - rg*i_LineD) / Lg;
    dx(12) = (v_gfm_Q - V_inf_Q - omega0_nominal*Lg*i_LineD - rg*i_LineQ) / Lg;
end

function x_lim = calc_limiter(x, lim)
    min_lim = lim(1);
    max_lim = lim(2);

    if x >= max_lim
        x_lim = max_lim;
    elseif x <= min_lim
        x_lim = min_lim;
    else
        x_lim = x;
    end
end