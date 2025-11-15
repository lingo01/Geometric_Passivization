function [RGrid, IGrid, wPlot] = func_PRPregion_EX2(R, I, sigma, nR, nI)
    % positive real region
    % Number of grid points (can be adjusted as needed)
    if nargin == 3
        nR = 30; 
        nI = 30;
    elseif nargin == 4
        nI = nR;
    end
    % Non-uniform sampling parameter: gamma>1 makes sampling denser near 0
    gamma = 2.5;

    % Save original range (to avoid overwriting input variable names)
    R_range = R; I_range = I;

    % Generate R axis (clustered at 0)
    uR = linspace(-1, 1, nR);
    R = zeros(size(uR));
    Rmin = R_range(1); Rmax = R_range(2);
    posR = uR >= 0;
    % Positive half-axis (0 to Rmax) uses u^gamma mapping, negative half-axis is symmetric to Rmin
    if Rmax > 0
        R(posR) = (uR(posR) .^ gamma) * Rmax;
    else
        R(posR) = uR(posR) * Rmax; % If upper bound is not positive, revert to linear mapping
    end
    if Rmin < 0
        R(~posR) = - (abs(uR(~posR)) .^ gamma) * abs(Rmin);
    else
        R(~posR) = uR(~posR) * Rmin; % If lower bound is not negative, revert to linear mapping
    end

    % Generate I axis (clustered at 0)
    uI = linspace(-1, 1, nI);
    I = zeros(size(uI));
    Imin = I_range(1); Imax = I_range(2);
    posI = uI >= 0;
    if Imax > 0
        I(posI) = (uI(posI) .^ gamma) * Imax;
    else
        I(posI) = uI(posI) * Imax;
    end
    if Imin < 0
        I(~posI) = - (abs(uI(~posI)) .^ gamma) * abs(Imin);
    else
        I(~posI) = uI(~posI) * Imin;
    end

    [RGrid, IGrid] = meshgrid(R, I);
    
    denom = sigma * (RGrid.^2 + IGrid.^2);
    
    wPlot = nan(size(denom));
    mask = (IGrid < 0) & (denom > 0);
    wBoundary = zeros(size(denom));
    wBoundary(mask) = -IGrid(mask) ./ denom(mask);
    
    wBoundary(wBoundary > 100) = 100;
    
    wPlot(mask) = wBoundary(mask);
end