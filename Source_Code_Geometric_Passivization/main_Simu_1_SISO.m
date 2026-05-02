clear all; clc;

%%
s = tf('s');

T = 0.02;
M = 0.3;
d = 0.5;
tau = 0.1;
k = 0.5;

w = 10.^[-5:0.001:5];


G_1 = 1 / (tau*s+k); 
G_2 = 1 / (s*(M*s+d));
G_3 = 1 / ((T*s+1) * (M*s+d));

alpha = 0.5;
tr = 1;
te = 10000000;

G_5 = alpha / s * (1/s - (tr+1/s)*(1-tr/2*s)/(1+tr/2*s) + tr*(1-tr/2*s)/(1+tr/2*s)...
    - tr*(1-te/2*s)/(1+te/2*s));

%%
[R_G1, I_G1] = nyquist(G_1);  R_G1 = squeeze(R_G1); I_G1 = squeeze(I_G1); 
[R_G2, I_G2] = nyquist(G_2);  R_G2 = squeeze(R_G2); I_G2 = squeeze(I_G2); 
[R_G3, I_G3] = nyquist(G_3);  R_G3 = squeeze(R_G3); I_G3 = squeeze(I_G3); 
[R_G5, I_G5] = nyquist(G_5);  R_G5 = squeeze(R_G5); I_G5 = squeeze(I_G5); 
[R_sG5, I_sG5] = nyquist(s*G_5);  R_sG5 = squeeze(R_sG5); I_sG5 = squeeze(I_sG5); 

% (a)
figure(1);
subplot(3,2,1); title('(a)'); hold on; box on;
plot(R_G1, I_G1, 'LineWidth', 1.2, 'Color', [217,083,025]/255);

theta = linspace(pi, 2*pi, 200);
xc = 0.5 + 0.5*cos(theta);
yc = 0   + 0.5*sin(theta);
patch(xc, yc, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

theta = linspace(2*pi, pi, 200);
xc = 1.5 + 1.5*cos(theta);
yc = 0   + 1.5*sin(theta);
patch(xc, yc, 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
xlim([0, 3]); ylim([-3, 0]); 
xlabel('$\Re(G)$','Interpreter','latex');
ylabel('$\Im(G)$','Interpreter','latex');

ax = gca;
ax.FontName = 'Times New Roman';
ax.XAxis.FontName = 'Times New Roman';
ax.YAxis.FontName = 'Times New Roman';
try
    ax.ZAxis.FontName = 'Times New Roman';
end
ax.TickLabelInterpreter = 'tex';

% (b)
subplot(3,2,2); title('(b)'); hold on; box on;

plot(R_G2, I_G2, 'LineWidth', 1.2, 'Color', [217,083,025]/255);
theta = linspace(pi, 2*pi, 200);
xc = 1.5 + 1.5*cos(theta);
yc = 0   + 1.5*sin(theta);
patch(xc, yc, 'b', 'FaceAlpha', 0.25, 'EdgeColor', 'none');
xlim([-3, 3]); ylim([-6, 0]); 
xlabel('$\Re(G)$','Interpreter','latex');
ylabel('$\Im(G)$','Interpreter','latex');

ax = gca;
ax.FontName = 'Times New Roman';
ax.XAxis.FontName = 'Times New Roman';
ax.YAxis.FontName = 'Times New Roman';
try
    ax.ZAxis.FontName = 'Times New Roman';
end
ax.TickLabelInterpreter = 'tex';

% (c)
subplot(3,2,3); title('(c)'); hold on; box on;
plot(R_G3, I_G3, 'LineWidth', 1.2, 'Color', [217,083,025]/255);

theta = linspace(pi, 2*pi, 200);
xc = 1.5 + 1.5*cos(theta);
yc = 0   + 1.5*sin(theta);
patch(xc, yc, 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
xlim([-0.5, 3]); ylim([-3, 0]); 
xlabel('$\Re(G)$','Interpreter','latex');
ylabel('$\Im(G)$','Interpreter','latex');

ax = gca;
ax.FontName = 'Times New Roman';
ax.XAxis.FontName = 'Times New Roman';
ax.YAxis.FontName = 'Times New Roman';
try
    ax.ZAxis.FontName = 'Times New Roman';
end
ax.TickLabelInterpreter = 'tex';

% (e)
subplot(3,2,5); 
title('(e)'); hold on; box on;
plot(R_G5, I_G5, 'LineWidth', 1.2, 'Color', [217,083,025]/255);
plot(R_G5, -I_G5, '--', 'LineWidth', 1.2, 'Color', [217,083,025]/255);

% (f)
subplot(3,2,6); 
title('(f)'); hold on; box on;
plot(R_sG5, I_sG5, 'LineWidth', 1.2, 'Color', [217,083,025]/255);
plot(R_sG5, -I_sG5, '--', 'LineWidth', 1.2, 'Color', [217,083,025]/255);

theta = linspace(0, 2*pi, 200);
a = 1/(2*1/(3*alpha*tr));
xc = a + a*cos(theta);
yc = 0 + a*sin(theta);
patch(xc, yc, 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
% xlim([0,1]); ylim([-1,1]/2)
xlabel('$\Re(G)$','Interpreter','latex');
ylabel('$\Im(G)$','Interpreter','latex');

ax = gca;
ax.FontName = 'Times New Roman';
ax.XAxis.FontName = 'Times New Roman';
ax.YAxis.FontName = 'Times New Roman';
try
    ax.ZAxis.FontName = 'Times New Roman';
end
ax.TickLabelInterpreter = 'tex';

%% (d)
sigma = [-0.5, -0.2, 0, 0.2, 0.5];

figure(1);
subplot(3,2,4); title('(d)'); hold on; box on;
nLines = length(sigma);
start = [0 0 0.4];
stop  = [0.6 0.85 1];
t = linspace(0,1,nLines)';
cmap = repmat(start, nLines, 1) + t*(stop - start);

legend_labels = cell(1, nLines);

for rr = 1:nLines

    Gtilde = G_3 / (1-sigma(rr)*G_3);

    [mag, phase] = bode(Gtilde, w);

    mag = squeeze(mag);
    phase = squeeze(phase);

    plot(log10(w), phase, 'Color', cmap(rr,:), 'LineWidth', 1.2);
    
    legend_labels{rr} = sprintf('\\sigma=%.1f', sigma(rr));
    
    phase_minus90_idx = find(phase <= -90, 1, 'first');
    if ~isempty(phase_minus90_idx)
        freq_minus90 = w(phase_minus90_idx);
        fprintf('sigma=%.1f, cross freq: %.4f rad/s (log10(w)=%.4f)\n',sigma(rr), freq_minus90, log10(freq_minus90));
    else
        fprintf('sigma=%.1f, cross freq: NaN\n', sigma(rr));
    end
end

legend(legend_labels);
ylim([-180, 0]);
xlabel('$\log (\omega)$','Interpreter','latex');
ylabel('$\arg (G)$','Interpreter','latex');

ax = gca;
ax.FontName = 'Times New Roman';
ax.XAxis.FontName = 'Times New Roman';
ax.YAxis.FontName = 'Times New Roman';
try
    ax.ZAxis.FontName = 'Times New Roman';
end
ax.TickLabelInterpreter = 'tex';