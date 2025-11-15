% main_Illu_2_Ctrl.m
% Draw a Nichols chart of a loop transfer function and overlay the shaded
% region defined by: for theta in [-pi/2, pi/2], allowed linear gain <= cos(theta)/sigma.

clear; close all; clc;

%% Settings
s = tf('s');

% Example loop transfer function L(s); replace with your own as needed
tau = 0.1;
k = 0.5;

L = 1 / (tau*s+k); 

% Output feedback passivity index sigma (change as needed)
sigma = 0.1;

% Frequency vector for Nichols
w = logspace(-5, 5, 1000);


%% Nichols plot
figure(6);
[mag, phase] = nichols(L, w); 
plot(squeeze(phase), squeeze(20*log10(mag)), 'LineWidth', 1.2, 'Color', [217,083,025]/255);
ngrid;
hold on;

%% Positive Real Part Region (Nichols-domain equivalent)

% Use degrees for Nichols phase axis
theta_deg = linspace(-90, 90, 2000);

% Compute linear gain limit = cos(theta)/sigma, only where cos>0
cosvals = cosd(theta_deg);
mask_pos = cosvals > 0;
limit_linear = nan(size(cosvals));
limit_linear(mask_pos) = cosvals(mask_pos) / sigma;
limit_db = 20*log10(limit_linear);

% Fill polygon down to a low dB floor for plotting
mag_min_db = -140;
limit_db_filled = limit_db;
limit_db_filled(~mask_pos) = mag_min_db;

x_patch = [theta_deg, fliplr(theta_deg)];
y_patch = [limit_db_filled, repmat(mag_min_db, size(theta_deg))];

% Draw filled patch using solid blue and specified transparency, no edge
hPatch = patch(x_patch, y_patch, 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% (no colored line on the edge per request)

%% Plot Settings

% Formatting
% title('(a)', 'FontName', 'Times New Roman', 'Interpreter', 'none');
xlabel(['$\angle G$',' /deg'],'Interpreter','latex');
ylabel(['$20\log|G|$',' /dB'],'Interpreter','latex');

% Set finite x-limits to avoid excessive phase tick calculation
xlim([-200, 0]);
ylim([-50, Inf]);

% Ensure all existing figure objects use Times New Roman where supported
allObjs = findall(gcf);
for k = 1:numel(allObjs)
	try
		set(allObjs(k), 'FontName', 'Times New Roman');
	catch
		% object does not support FontName - ignore
	end
	% Also set Interpreter for text objects if they use LaTeX and you want
	% to keep LaTeX, we don't change Interpreter here to avoid breaking labels
end

hold off;
% End of file
