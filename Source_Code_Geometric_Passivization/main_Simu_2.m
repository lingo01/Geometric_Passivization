clear all; clc;

%%
s = tf('s');

T = 0.02;
M = 0.3;
d = 0.5;
tau = 0.1;
k = 0.5;

w = 10.^[-4:0.01:5];


G_1 = 1 / (tau*s+k); 
G_2 = 1 / (s*(M*s+d));
G_3 = 1 / ((T*s+1) * (M*s+d));


%%
[R_G2, I_G2] = nyquist(G_2);  R_G2 = squeeze(R_G2); I_G2 = squeeze(I_G2); 
[R_G3, I_G3] = nyquist(G_3);  R_G3 = squeeze(R_G3); I_G3 = squeeze(I_G3); 

% (a)
figure(2);
subplot(1,2,1); title('(a)'); hold on; box on;
plot(R_G2, I_G2, 'LineWidth', 1.2, 'Color', [217,083,025]/255);
xlim([-3, 3]); ylim([-6, 0.2]); 
xlabel('$\Re(G)$','Interpreter','latex');
ylabel('$\Im(G)$','Interpreter','latex');

xl = xlim; yl = ylim; 
patch([xl(1), xl(2), xl(2), xl(1)], [yl(1), yl(1), 0, 0], 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');

ax = gca;
ax.FontName = 'Times New Roman';
ax.XAxis.FontName = 'Times New Roman';
ax.YAxis.FontName = 'Times New Roman';
try
    ax.ZAxis.FontName = 'Times New Roman';
end
ax.TickLabelInterpreter = 'tex';

% (b)
subplot(1,2,2); title('(b)'); hold on; box on;
plot(R_G3, I_G3, 'LineWidth', 1.2, 'Color', [217,083,025]/255);
xlim([-0.5, 3]); ylim([-3, 0.1]); 
xlabel('$\Re(G)$','Interpreter','latex');
ylabel('$\Im(G)$','Interpreter','latex');

xl = xlim; yl = ylim; 
patch([xl(1), xl(2), xl(2), xl(1)], [yl(1), yl(1), 0, 0], 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');

ax = gca;
ax.FontName = 'Times New Roman';
ax.XAxis.FontName = 'Times New Roman';
ax.YAxis.FontName = 'Times New Roman';
try
    ax.ZAxis.FontName = 'Times New Roman';
end
ax.TickLabelInterpreter = 'tex';

%% Save figure(2) as a vector PDF
% This tries to use exportgraphics (R2020a+) to ensure true vector output.
% If exportgraphics is not available, fallback to print with the painters renderer.
% try
%     fig = figure(2);
%     % Ensure the figure uses a vector-friendly renderer
%     set(fig, 'Renderer', 'painters');
%     outputFile = fullfile(pwd, 'Simulation2.pdf');
% 
%     if exist('exportgraphics', 'file') == 2
%         % exportgraphics supports ContentType='vector'
%         exportgraphics(fig, outputFile, 'ContentType', 'vector', 'BackgroundColor', 'none');
%     else
%         % Fallback for older MATLAB versions
%         print(fig, outputFile, '-painters', '-dpdf', '-bestfit');
%     end
% 
%     fprintf('Saved figure(1) to %s\n', outputFile);
% catch ME
%     warning('Failed to save figure(1) as PDF: %s', ME.message);
% end