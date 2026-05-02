function [kp, ki, feasible] = func_solve_PI_passivity(sys, sigma, wc, param_a, param_b)
    % Extract the minimal state-space realization of the plant
    [Ap, Bp, Cp, Dp] = ssdata(sys);
    np = size(Ap, 1);
    
    % Construct the full-dimensional augmented constant matrix when PI controller is in observable canonical form
    % AL dimension: (np+1) x (np+1)
    AL = [-param_b/param_a, zeros(1, np); Bp, Ap];
    % CL dimension: 1 x (np+1)
    CL = [Dp, Cp];
    
    % Construct global projection matrices MB and MD for the affine parameters
    % so that BL(rho) = MB * rho, DL(rho) = MD * rho
    MB = [0, 1/param_a; Bp, zeros(np, 1)];
    MD = [Dp, 0];
    
    nL = np + 1; % Augmented system dimension

    % Frequency scaling (reduces the condition number for interior-point solvers)
    scaling_factor = 1e3;
    AL = AL / scaling_factor;
    CL = CL / sqrt(scaling_factor);
    MB = MB / sqrt(scaling_factor);
    MD = MD;
    wc = wc / scaling_factor;
    
    % Declare YALMIP decision variables
    P = sdpvar(nL, nL, 'symmetric');
    Q = sdpvar(nL, nL, 'symmetric');
    rho = sdpvar(2, 1); % Parameter vector: rho = [kp; ki]
    
    % Construct local block matrices for LMI
    Xi11 = AL*P + P*AL' - AL*Q*AL' + wc^2*Q;
    Xi12 = P*CL' - AL*Q*CL' - (1/(2*sigma)) * MB * rho;
    Xi13 = MB * rho;
    Xi22 = -CL*Q*CL' - (1/sigma) * MD * rho;
    Xi23 = MD * rho;
    
    % Assemble global matrix inequality
    Xi = [Xi11, Xi12, Xi13;
          Xi12', Xi22, Xi23;
          Xi13', Xi23', -1];
          
    % Set tolerance to prevent numerical singularity
    tol = 1e-6; 
    
    % Define constraints: LMI negative definite, Q positive definite, PID parameters non-negative
    Constraints = [Xi <= -tol*eye(size(Xi)), Q >= tol*eye(nL), rho >= 0];
    
    % Configure solver parameters (turn off verbose output)
    ops = sdpsettings('solver', 'sedumi', 'verbose', 0); 
    
    % Solve (no objective function means finding a feasible solution)
    sol = optimize(Constraints, 0, ops);
    
    if sol.problem == 0
        feasible = true;
        kp = double(rho(1));
        ki = double(rho(2));
        fprintf('LMI solved: kp = %.4f, ki = %.4f\n', kp, ki);
    else
        feasible = false;
        kp = NaN;
        ki = NaN;
        fprintf('LMI not solved.\n');
    end
end
