function results = func_checkNyquist(Gs, ReVals, ImVals, omega, sigma, tol, figHandle)
%FUNC_CHECKNYQUISTPR Check Nyquist curves against PR-region: w*sigma*(Re^2+Im^2)+Im <= 0
%   results = func_checkNyquistPR(Gs, ReVals, ImVals, omega, sigma, tol, figHandle)
%
%   Inputs:
%     Gs       - cell array of LTI models (for local refinement), or empty {}
%     ReVals   - cell array of real parts for each system
%     ImVals   - cell array of imag parts for each system
%     omega    - frequency vector (same length for each Re/Im)
%     sigma    - scalar sigma used in inequality
%     tol      - numeric tolerance (optional, default 1e-9)
%     figHandle- figure handle to plot violation markers (optional)
%
%   Output:
%     results  - struct array (length = number of systems checked) with fields:
%       inside    - logical (1 or true means ALL sampled points satisfy the PR inequality
%                   w*sigma*(Re^2+Im^2) + Im <= 0; 0/false means there exists at least one
%                   sampled point that violates the inequality)
%       maxv      - the maximum value of v = w*sigma*(Re^2+Im^2) + Im over the sampled points
%                   (useful to see how far from the boundary the curve is; negative means
%                   comfortably inside, small positive may be numerical)
%       violIdx   - indices (into the per-system omega/Re/Im vectors) where v > tol
%       violOmega - frequencies corresponding to violIdx (empty if no violations)
%       localMax  - if a violation exists and a corresponding LTI model was provided in Gs,
%                   this is the maximum v found after a local fine re-sampling around the
%                   first violation (NaN if not computed)

%   Note on interpretation:
%     - inside == 1 (true): all sampled points passed the test (<= tol). This means the
%       Nyquist curve, at the sampling resolution used, lies inside the positive-real
%       region. It does NOT strictly prove the continuous curve is entirely inside, but
%       is a practical numeric check. For higher assurance, increase sampling density or
%       rely on the local refinement results.
%     - inside == 0 (false): there exists at least one sampled point that violates the
%       inequality (v > tol). The function will report the violating frequencies and
%       perform one local refinement near the first violating point if the LTI model is
%       provided in Gs.


if nargin < 6 || isempty(tol)
    tol = 1e-9;
end
if nargin < 7
    figHandle = [];
end

n = numel(ReVals);
results = struct('inside', cell(1,n), 'maxv', [], 'violIdx', [], 'violOmega', [], 'localMax', []);

% colors for markers (match main script)
colors = {[217,83,25]/255, [128,0,32]/255, [0,0,139]/255};

checkValue = @(Re,Im,w,sigma) w .* sigma .* (Re.^2 + Im.^2) + Im;

for i = 1:n
    Re = ReVals{i}; Im = ImVals{i};
    % determine omega for this system: support numeric vector (shared) or cell array per-system
    if iscell(omega)
        if numel(omega) >= i && ~isempty(omega{i})
            omega_i = omega{i};
        else
            % fallback: create a linspace in [0,100] with same length as Re
            omega_i = linspace(0, 100, numel(Re));
        end
    else
        omega_common = omega;
        if numel(omega_common) == numel(Re)
            omega_i = omega_common;
        else
            % lengths differ: create omega_i by mapping omega_common range to match Re length
            if ~isempty(omega_common)
                omega_i = linspace(min(omega_common), max(omega_common), numel(Re));
            else
                omega_i = linspace(0, 100, numel(Re));
            end
        end
    end

    v = checkValue(Re, Im, omega_i, sigma);
    results(i).maxv = max(v);
    results(i).violIdx = find(v > tol);
    if isempty(results(i).violIdx)
        results(i).inside = true;
        results(i).violOmega = [];
        results(i).localMax = NaN;
    else
        results(i).inside = false;
        % safe assignment of violOmega using per-system omega_i
        results(i).violOmega = omega_i(results(i).violIdx);
        % plotting violations if figure provided
        if ~isempty(figHandle) && ishandle(figHandle)
            figure(figHandle); hold on;
            scatter3(Re(results(i).violIdx), Im(results(i).violIdx), omega_i(results(i).violIdx), 60, colors{min(i,numel(colors))}, 'filled');
            hold off;
        end
        % local refinement around first violation if Gs provided
        results(i).localMax = NaN;
        if ~isempty(Gs) && numel(Gs) >= i && ~isempty(Gs{i})
            Gsys = Gs{i};
            idx = results(i).violIdx(1);
            lo = max(idx-1, 1); hi = min(idx+1, numel(omega_i));
            omega_fine = linspace(omega_i(lo), omega_i(hi), 501);
            H_fine = squeeze(freqresp(Gsys, omega_fine));
            v_fine = checkValue(real(H_fine), imag(H_fine), omega_fine, sigma);
            results(i).localMax = max(v_fine);
        end
    end
end

end
