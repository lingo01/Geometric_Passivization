clear all; clc; close all;

%%
s = tf('s');

T = 0.02;
M = 0.3;
d = 0.5;
tau = 0.1;
k = 0.5;


G_1 = 1 / (tau*s+k);
G_2 = 1 / (s*(M*s+d));
G_3 = 1 / ((T*s+1) * (M*s+d));

% parameters for EX3 region
epsilon = 0.05; % user-specified constant (change as needed)
sigma = 0.1; % sigma parameter (shared with EX2-style region)

NGrid = 200;
NFreq = 1000;

%% EX3 region

[RGrid, IGrid, wPlot] = func_PRPregion_EX3([-15, 15], [-15, 0.1], epsilon, sigma, 2*NGrid);

figure(5); subplot(1,2,1);
surf(RGrid, IGrid, wPlot, 'FaceColor', 'b', 'FaceAlpha', 0.35, 'EdgeColor', 'none');
hold on;

w0 = zeros(size(RGrid));
w0(~(IGrid <= 0)) = NaN;
surf(RGrid, IGrid, w0, 'FaceColor', 'b', 'FaceAlpha', 0.08, 'EdgeColor', 'none');

xlabel('$\Re(G(j\omega))$','Interpreter','latex');
ylabel('$\Im(G(j\omega))$','Interpreter','latex');
zlabel('$\omega$','Interpreter','latex');
xlim([-15 15]); ylim([-15 0.1]); zlim([0 80]);
view(3); grid on; box on;
camlight headlight; lighting gouraud;
% Ensure numeric tick labels use Times New Roman (explicit axis properties)
ax = gca;
ax.FontName = 'Times New Roman';
ax.XAxis.FontName = 'Times New Roman';
ax.YAxis.FontName = 'Times New Roman';
try
    ax.ZAxis.FontName = 'Times New Roman';
end
ax.TickLabelInterpreter = 'tex';

%% Nyquist plot of transfer functions (overlay)

[G1val, ReG1val, ImG1val, omega_nyq] = func_Nyquist(G_1, NFreq);

h1 = plot3(ReG1val, ImG1val, omega_nyq, 'LineWidth', 2, 'Color', [217, 83, 25]/255);
plot3(ReG1val(1), ImG1val(1), omega_nyq(1), 'o', 'MarkerFaceColor', [217, 83, 25]/255, 'MarkerEdgeColor', 'none');

[G2val, ReG2val, ImG2val, omega_nyq] = func_Nyquist(G_2, NFreq);

h2 = plot3(ReG2val, ImG2val, omega_nyq, 'LineWidth', 2, 'Color', [128, 0, 32]/255);
plot3(ReG2val(1), ImG2val(1), omega_nyq(1), 'o', 'MarkerFaceColor', [128, 0, 32]/255, 'MarkerEdgeColor', 'none');

[G3val, ReG3val, ImG3val, omega_nyq] = func_Nyquist(G_3, NFreq);

h3 = plot3(ReG3val, ImG3val, omega_nyq, 'LineWidth', 2, 'Color', [0, 0, 139]/255);
plot3(ReG3val(1), ImG3val(1), omega_nyq(1), 'o', 'MarkerFaceColor', [0, 0, 139]/255, 'MarkerEdgeColor', 'none');

% legends
legend([h1, h2, h3], {'$G_1$', '$G_2$', '$G_3$'}, 'Location', 'best','Interpreter','latex', 'FontName', 'Times New Roman');
title('(a)', 'FontName', 'Times New Roman')
hold off;

%% comparison: EX2 vs EX3 for same sigma
figure(5); subplot(1,2,2);
hold on;

% EX2 region (assumes func_PRregion exists and matches EX2 behavior)
[R2, I2, w2] = func_PRPregion_EX2([-15, 15], [-15, 0.1], sigma, NGrid);
    color2 = [0.0,  0.45, 0.74];
C2 = repmat(reshape(color2, [1 1 3]), size(w2,1), size(w2,2));
hEx2 = surf(R2, I2, w2, C2, 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 0.5);

% EX3 region (same sigma)
[R3, I3, w3] = func_PRPregion_EX3([-15, 15], [-15, 0.1], epsilon, sigma, NGrid);
    color3 = [0.85, 0.85, 0.85];
C3 = repmat(reshape(color3, [1 1 3]), size(w3,1), size(w3,2));
hEx3 = surf(R3, I3, w3, C3, 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 0.5);

xlabel('$\Re(G(j\omega))$','Interpreter','latex');
ylabel('$\Im(G(j\omega))$','Interpreter','latex');
zlabel('$\omega$','Interpreter','latex');
xlim([-15 15]); ylim([-15 0.1]); zlim([0 80]);
view(3); box on;
camlight headlight; lighting none;
% Ensure numeric tick labels use Times New Roman (explicit axis properties)
ax = gca;
ax.FontName = 'Times New Roman';
ax.XAxis.FontName = 'Times New Roman';
ax.YAxis.FontName = 'Times New Roman';
try
    ax.ZAxis.FontName = 'Times New Roman';
end
ax.TickLabelInterpreter = 'tex';

legend([hEx2, hEx3], {'Ex2: $\sigma=0.1$','Ex3: ($\sigma,\varepsilon)=(0.1,0.02)$'}, 'Interpreter', 'latex', 'Location', 'eastoutside', 'FontSize', 10, 'FontName', 'Times New Roman');
title('(b)', 'Interpreter', 'latex', 'FontName', 'Times New Roman')

set(gcf,'Renderer','opengl');

hold off;

% Save figure(5) as a high-quality PNG
try
    fig = figure(5);
    % Prefer a modern renderer for high-quality raster output
    set(fig, 'Renderer', 'opengl');

    outputFile = fullfile(pwd, 'Simulation3_EX3_easy.png');
    dpi = 600; % change to 300 if file size is too large

    if exist('exportgraphics', 'file') == 2
        % exportgraphics supports a Resolution option
        exportgraphics(fig, outputFile, 'Resolution', dpi);
        fprintf('Saved high-quality PNG to %s (dpi=%d) using exportgraphics.\n', outputFile, dpi);
    else
        % Fallback for older MATLAB versions
        print(fig, outputFile, '-dpng', sprintf('-r%d', dpi));
        fprintf('Saved high-quality PNG to %s (dpi=%d) using print.\n', outputFile, dpi);
    end
catch ME
    warning('Failed to save PNG: %s', ME.message);
end