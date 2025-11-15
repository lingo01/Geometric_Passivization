clear all; clc;

%%
s = tf('s');

T = 0.02;
M = 0.3;
d = 0.5;
tau = 0.1;
k = 0.5;
C = 0.1;


G_1 = 1 / (tau*s+k); 
G_2 = 1 / (s*(M*s+d));
G_3 = 1 / ((T*s+1) * (M*s+d));
G_c = C / (T*s+1);
G_4 = [G_3, G_c; G_c, G_1];


%% Rayleigh quotient of G
w_seq_note = [0.1, 0.5, 1, 10, 20, 50, 100];

w_seq = 10.^[-3:0.1:2];

figure(4); clf; hold on; axis equal;

% Use randomized sampling + convex hull to produce a filled, simply-connected
% numerical range for each frequency (avoids self-intersecting patch issues).
nFreq = length(w_seq);
Nsample = 8000; % samples per frequency (adjust for speed/quality)

% stronger perceptual gradient: create a custom red -> purple -> blue colormap
tcol = linspace(0,1,nFreq)';
c1 = [1, 0, 0];      % red
c2 = [0.6, 0, 0.6];  % purple
c3 = [0, 0, 1];      % blue
cmap = zeros(nFreq,3);
for ii = 1:nFreq
    tt = tcol(ii);
    if tt <= 0.5
        s = tt / 0.5;
        cmap(ii,:) = (1-s)*c1 + s*c2;
    else
        s = (tt-0.5) / 0.5;
        cmap(ii,:) = (1-s)*c2 + s*c3;
    end
end

% pick target index closest to w = 1e-1 for dashed-circle highlighting
[~, idx_target] = min(abs(w_seq - 10^0.7));

for rr = 1:nFreq
    w = w_seq(rr);
    % evaluate transfer-function matrix G_4 at s = j*w (should be 2x2 numeric)
    G4w = evalfr(G_4, 1j*w);

    % Random sampling (vectorized): generate complex random vectors
    X = (randn(2, Nsample) + 1j * randn(2, Nsample)); % 2 x N
    Y = G4w * X;                                    % 2 x N
    num = sum(conj(X) .* Y, 1);                    % 1 x N
    den = sum(conj(X) .* X, 1);                    % 1 x N
    qs = num ./ den;                                % 1 x N Rayleigh quotients

    % Convex hull of sampled points to get a clean boundary
    kidx = convhull(real(qs), imag(qs));

    % choose color from colormap and alpha (more pronounced gradient)
    color = cmap(rr, :);
    frac = (rr-1) / max(1, (nFreq-1));
    alphaVal = 0.08 + 0.7 * frac;     % transparency increases with freq index

    % fill the convex hull; conditional edge styling
    if rr == idx_target
        % target: black dashed edge to highlight
        patch(real(qs(kidx)), imag(qs(kidx)), color, 'FaceAlpha', alphaVal, ...
            'EdgeColor', 'k', 'LineWidth', 2);
    else
        % non-target: faint/no edge
        patch(real(qs(kidx)), imag(qs(kidx)), color, 'FaceAlpha', alphaVal, ...
            'EdgeColor', 'none');
    end
end

%% positive real region
t = linspace(0, 2*pi, 360);
circ = 1.5 + 1.5 * exp(1j * t);
plot(real(circ), imag(circ), '--', 'LineWidth', 2, 'Color', [217,083,025]/255);

xlabel('$\Re (\rho_G)$','Interpreter','latex');
ylabel('$\Im (\rho_G)$','Interpreter','latex');
xlim([-0.2, 3]); ylim([-1.5, 0]);
box on; grid off;

ax = gca;
ax.FontName = 'Times New Roman';
ax.XAxis.FontName = 'Times New Roman';
ax.YAxis.FontName = 'Times New Roman';
try
    ax.ZAxis.FontName = 'Times New Roman';
end
ax.TickLabelInterpreter = 'tex';

hold off;


%% Save figure(4) as a vector PDF
% This tries to use exportgraphics (R2020a+) to ensure true vector output.
% If exportgraphics is not available, fallback to print with the painters renderer.
% try
%     fig = figure(4);
%     % Ensure the figure uses a vector-friendly renderer
%     set(fig, 'Renderer', 'painters');
%     outputFile = fullfile(pwd, 'Simulation1_MIMO.pdf');
% 
%     if exist('exportgraphics', 'file') == 2
%         % exportgraphics supports ContentType='vector'
%         exportgraphics(fig, outputFile, 'ContentType', 'vector', 'BackgroundColor', 'none');
%     else
%         % Fallback for older MATLAB versions
%         print(fig, outputFile, '-painters', '-dpdf', '-bestfit');
%     end
% 
%     fprintf('Saved figure(4) to %s\n', outputFile);
% catch ME
%     warning('Failed to save figure(1) as PDF: %s', ME.message);
% end

