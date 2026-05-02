function func_plot_PI_feasible_region(sys, sigma, wc, param_a, param_b)
    % Extract matrices and augmented construction (consistent with solver function)
    [Ap, Bp, Cp, Dp] = ssdata(sys);
    np = size(Ap, 1);
    AL = [-param_b/param_a, zeros(1, np); Bp, Ap];
    CL = [Dp, Cp];
    MB = [0, 1/param_a; Bp, zeros(np, 1)];
    MD = [Dp, 0];
    nL = np + 1;

    % Frequency scaling (reduces the condition number for interior-point solvers)
    scaling_factor = 1e3;
    AL = AL / scaling_factor;
    CL = CL / sqrt(scaling_factor);
    MB = MB / sqrt(scaling_factor);
    MD = MD;
    wc = wc / scaling_factor;
    
    P = sdpvar(nL, nL, 'symmetric');
    Q = sdpvar(nL, nL, 'symmetric');
    rho = sdpvar(2, 1);
    
    Xi11 = AL*P + P*AL' - AL*Q*AL' + wc^2*Q;
    Xi12 = P*CL' - AL*Q*CL' - (1/(2*sigma)) * MB * rho;
    Xi13 = MB * rho;
    Xi22 = -CL*Q*CL' - (1/sigma) * MD * rho;
    Xi23 = MD * rho;
    
    Xi = [Xi11, Xi12, Xi13; Xi12', Xi22, Xi23; Xi13', Xi23', -1];
    tol = 1e-6;
    
    % Add upper bound on parameters to prevent unbound feasible region in some directions
    param_max = 50; 
    Constraints = [Xi <= -tol*eye(size(Xi)), Q >= tol*eye(nL),...
                   rho >= 0, rho <= param_max];
               
    ops = sdpsettings('solver', 'sedumi', 'verbose', 0);
    
    % Emit optimization vectors in various angles to obtain convex hull boundary
    angles = linspace(0, 2*pi, 40);
    kp_boundary = [];
    ki_boundary = [];
    
    for theta = angles
        % Set objective to maximize in the theta direction
        obj = -(cos(theta)*rho(1) + sin(theta)*rho(2)); 
        sol = optimize(Constraints, obj, ops);
        
        if sol.problem == 0 || sol.problem == 4 % 0 is optimal, 4 is unbounded (truncated by param_max)
            kp_boundary(end+1) = double(rho(1));
            ki_boundary(end+1) = double(rho(2));
        end
    end
    
    % Plot convex hull region
    if ~isempty(kp_boundary)
        figure;
        K = convhull(kp_boundary, ki_boundary);
        fill(kp_boundary(K), ki_boundary(K), [0.2 0.6 0.8], 'FaceAlpha', 0.4, 'EdgeColor', 'b', 'LineWidth', 1.5);
        xlabel('$k_p$', 'Interpreter', 'latex');
        ylabel('$k_i$', 'Interpreter', 'latex');
        title(sprintf('PI parameter feasible region (\\sigma = %.1f, \\omega_c = %.1f)', sigma, wc));
        grid on;
    else
        disp('PI parameter feasible region does not exist.');
    end
end
