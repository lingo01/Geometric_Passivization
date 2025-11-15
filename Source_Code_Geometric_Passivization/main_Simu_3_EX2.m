clear all; clc;

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

sigma = 0.1;

NGrid = 200;
NFreq = 1000;

%% positive real region

[RGrid, IGrid, wPlot] = func_PRPregion_EX2([-2, 2], [-2, 0.1], sigma, 2*NGrid);

figure(3); subplot(1,2,1);
surf(RGrid, IGrid, wPlot, 'FaceColor', 'b', 'FaceAlpha', 0.35, 'EdgeColor', 'none');
hold on;

w0 = zeros(size(RGrid));
w0(~(IGrid <= 0)) = NaN;
surf(RGrid, IGrid, w0, 'FaceColor', 'b', 'FaceAlpha', 0.08, 'EdgeColor', 'none');

xlabel('$\Re(G(j\omega))$','Interpreter','latex');
ylabel('$\Im(G(j\omega))$','Interpreter','latex');
zlabel('$\omega$','Interpreter','latex');
xlim([-2 2]); ylim([-2 0.1]); zlim([0 100]);
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

%% Nyquist plot of transfer functions

% nyquist plot of G1

[G1val, ReG1val, ImG1val, omega_nyq] = func_Nyquist(G_1, NFreq);

h1 = plot3(ReG1val, ImG1val, omega_nyq, 'LineWidth', 2, 'Color', [217, 83, 25]/255);
plot3(ReG1val(1), ImG1val(1), omega_nyq(1), 'o', 'MarkerFaceColor', [217, 83, 25]/255, 'MarkerEdgeColor', 'none');

% nyquist plot of G_2

[G2val, ReG2val, ImG2val, omega_nyq] = func_Nyquist(G_2, NFreq);

h2 = plot3(ReG2val, ImG2val, omega_nyq, 'LineWidth', 2, 'Color', [128, 0, 32]/255);
plot3(ReG2val(1), ImG2val(1), omega_nyq(1), 'o', 'MarkerFaceColor', [128, 0, 32]/255, 'MarkerEdgeColor', 'none');


% nyquist plot of G_3

[G3val, ReG3val, ImG3val, omega_nyq] = func_Nyquist(G_3, NFreq);

h3 = plot3(ReG3val, ImG3val, omega_nyq, 'LineWidth', 2, 'Color', [0, 0, 139]/255);
plot3(ReG3val(1), ImG3val(1), omega_nyq(1), 'o', 'MarkerFaceColor', [0, 0, 139]/255, 'MarkerEdgeColor', 'none');

% legends
legend([h1, h2, h3], {'$G_1$', '$G_2$', '$G_3$'}, 'Location', 'best','Interpreter','latex', 'FontName', 'Times New Roman');
title('(a)', 'FontName', 'Times New Roman')
hold off;


%% check if Nyquist plot is within the positive real region
tol = 1e-9;
Gs = {G_1, G_2, G_3};
ReVals = {ReG1val, ReG2val, ReG3val};
ImVals = {ImG1val, ImG2val, ImG3val};
figHandle = figure(3);
results = func_checkNyquist(Gs, ReVals, ImVals, omega_nyq, sigma, tol, figHandle);


for i = 1:numel(results)
    fprintf('G%d: inside=%d, max v=%.3e, local refine max=%.3e\n', i, results(i).inside, results(i).maxv, results(i).localMax);
end

%% comparision of different passivity indices

sigma_seq = [0.1, 0.3, 1.0];

baseColors = [
    0.85, 0.85, 0.85;  % light gray
    % 0.6,  0.8,  1.0;   % light blue
    0.0,  0.45, 0.74;  % deeper blue
    0.49, 0.18, 0.56;  % purple
];

figure(3); subplot(1,2,2);
hold on;
hPR = gobjects(length(sigma_seq),1);
for rr = 1:length(sigma_seq)
    [RGrid, IGrid, wPlot] = func_PRPregion_EX2([-1, 1], [-1, 0.1], sigma_seq(rr), NGrid);

    color = baseColors(rr, :);
    alphaVal = 0.5;
    Cconst = repmat(reshape(color, [1 1 3]), size(wPlot,1), size(wPlot,2));
    h = surf(RGrid, IGrid, wPlot, Cconst, 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', alphaVal);
    hPR(rr) = h;
end

h3 = plot3(ReG3val, ImG3val, omega_nyq, 'LineWidth', 2, 'Color', [0, 0, 139]/255);
plot3(ReG3val(1), ImG3val(1), omega_nyq(1), 'o', 'MarkerFaceColor', [0, 0, 139]/255, 'MarkerEdgeColor', 'none');


xlabel('$\Re(G(j\omega))$','Interpreter','latex');
ylabel('$\Im(G(j\omega))$','Interpreter','latex');
zlabel('$\omega$','Interpreter','latex');
xlim([-1 1]); ylim([-1 0.1]); zlim([0 20]);
view(3); grid on; box on;
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

labels = arrayfun(@(s) sprintf('$\\sigma=%.3g$', s), sigma_seq, 'UniformOutput', false);
legend(hPR, labels, 'Interpreter', 'latex', 'Location', 'eastoutside', 'FontSize', 10, 'FontName', 'Times New Roman');
title('(b)', 'FontName', 'Times New Roman')

set(gcf,'Renderer','opengl');


hold off;




%% Save figure(3) as a vector PDF
% try
%     fig = figure(3);
%     % Ensure the figure uses a vector-friendly renderer
%     set(fig, 'Renderer', 'painters');
%     outputFile = fullfile(pwd, 'Simulation3.pdf');
% 
%     if exist('exportgraphics', 'file') == 2
%         % exportgraphics supports ContentType='vector'
%         exportgraphics(fig, outputFile, 'ContentType', 'vector', 'BackgroundColor', 'none');
%     else
%         % Fallback for older MATLAB versions
%         print(fig, outputFile, '-painters', '-dpdf', '-bestfit');
%     end
% 
%     fprintf('Saved figure(3) to %s\n', outputFile);
% end

% Save figure(3) as a high-quality PNG
try
    fig = figure(3);
    % Prefer a modern renderer for high-quality raster output
    set(fig, 'Renderer', 'opengl');

    outputFile = fullfile(pwd, 'Simulation3_EX2_easy.png');
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