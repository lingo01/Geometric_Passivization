function [RGrid, IGrid, wPlot] = func_PRPregion_EX3(R, I, e, sigma, nR, nI)
    % EX3 positive-region function
    % Inequality (updated): e*w*(R^2 - R) + (e*w*sigma*I^2 + I) <= 0
    % Solve for w (w >= 0): e*(R^2 - R + sigma*I^2)*w + I <= 0
    % If denom>0 and I<0 -> w <= -I/denom. If denom<=0 and I<=0 -> holds for all w>=0.

    if nargin == 4
        nR = 30;
        nI = 30;
    elseif nargin == 5
        nI = nR;
    end

    % Non-uniform sampling parameter (concentrate near 0)
    gamma = 2.5;

    % preserve original ranges
    R_range = R; I_range = I;

    % build R axis
    uR = linspace(-1, 1, nR);
    R = zeros(size(uR));
    Rmin = R_range(1); Rmax = R_range(2);
    posR = uR >= 0;
    if Rmax > 0
        R(posR) = (uR(posR) .^ gamma) * Rmax;
    else
        R(posR) = uR(posR) * Rmax;
    end
    if Rmin < 0
        R(~posR) = - (abs(uR(~posR)) .^ gamma) * abs(Rmin);
    else
        R(~posR) = uR(~posR) * Rmin;
    end

    % build I axis
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

    % For the requested inequality:
    % e*w*(sigma*R^2 - R) + (e*w*sigma*I^2 + I) <= 0
    % => e*(sigma*(R^2 + I^2) - R) * w + I <= 0
    % Let coeff = e*(sigma*(R^2 + I^2) - R)
    % If coeff > 0 and I <= 0 then w <= -I/coeff (finite upper bound).
    % If coeff <= 0 and I <= 0 then the inequality holds for all w>=0
    % (we leave those as NaN so the surface does not form an artificial cap).

    coeff = e * (sigma * (RGrid.^2 + IGrid.^2) - RGrid);
    wPlot = nan(size(coeff));

    % Case: finite upper bound when coeff > 0 and I < = 0
    maskUpper = (coeff > 0) & (IGrid <= 0);
    wBoundary = nan(size(coeff));
    wBoundary(maskUpper) = - IGrid(maskUpper) ./ coeff(maskUpper);
    % Only keep non-negative bounds (numerical safety)
    % wBoundary(wBoundary < 0) = NaN;
    
    wPlot(maskUpper) = wBoundary(maskUpper);

    % For coeff <= 0 and I <= 0 the inequality holds for all w>=0: keep NaN
    % For coeff > 0 and I > 0 there is no w>=0 that satisfies: keep NaN
end
