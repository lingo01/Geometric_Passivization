function [opt_sigma, opt_kp, opt_ki] = func_maximize_PI_passivity(sys, wc, param_a, param_b)
    % Set the search boundaries for the bisection method (physically, sigma won't be infinite)
    sigma_lower = 0.001; 
    sigma_upper = 100.0; 
    tol = 0.01; % Precision tolerance
    
    opt_sigma = NaN;
    opt_kp = NaN;
    opt_ki = NaN;
    
    fprintf('Starting bisection search for maximum passivity margin sigma...\n');
    
    % Check if the lower bound is feasible. If even a tiny sigma is infeasible, exit directly.
    [kp, ki, feasible] = func_solve_PI_passivity(sys, sigma_lower, wc, param_a, param_b);
    if ~feasible
        disp('Error: No feasible solution even for the minimum sigma.');
        return;
    end
    
    % Main bisection loop
    while (sigma_upper - sigma_lower) > tol
        sigma_mid = (sigma_lower + sigma_upper) / 2;
        fprintf('Testing sigma = %.4f... ', sigma_mid);
        
        % Call the single LMI solving function
        [kp, ki, feasible] = func_solve_PI_passivity(sys, sigma_mid, wc, param_a, param_b);
        
        if feasible
            % If the current sigma is feasible, the true maximum lies in the right interval
            fprintf('Feasible!\n');
            sigma_lower = sigma_mid; 
            opt_sigma = sigma_mid; % Record the current best parameters
            opt_kp = kp;
            opt_ki = ki;
        else
            % If the current sigma is infeasible, the true maximum lies in the left interval
            fprintf('Infeasible.\n');
            sigma_upper = sigma_mid;
        end
    end
    
    fprintf('=======================================\n');
    fprintf('Optimization complete!\n');
    fprintf('Maximum achievable passivity margin max_sigma = %.4f\n', opt_sigma);
    fprintf('Corresponding optimal PI parameters: kp = %.4f, ki = %.4f\n', opt_kp, opt_ki);
    fprintf('=======================================\n');
end