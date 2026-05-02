clear all; clc; close all;

%% settings
s = tf('s');

% transient simulation parameters
te = 10;
dt = 50e-6;
tspan = 0:dt:te;
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);

% Nyquist sampling frequency
w_nyq = logspace(-3, 5, 300);

% system parameters
[System_params, GFMInv_params] = setting_params(0);

% System strength settings
% stiff grid strength:  scaling_factor = 0.7
% normal grid strength: scaling_factor = 1.0
scaling_factor = 0.8; 
System_params.Lg = scaling_factor * System_params.Lg;
System_params.Xg = scaling_factor * System_params.Xg;

% compute equilibrium
x_equi = calc_equilibrium(GFMInv_params, System_params);

% update references with equilibrium values
v_od_eq = x_equi(1); v_oq_eq = x_equi(2); theta_gfm_eq = x_equi(9);
i_LineD_eq = x_equi(11); i_LineQ_eq = x_equi(12);
i_od_eq = i_LineD_eq * cos(theta_gfm_eq) + i_LineQ_eq * sin(theta_gfm_eq);
i_oq_eq = -i_LineD_eq * sin(theta_gfm_eq) + i_LineQ_eq * cos(theta_gfm_eq);
GFMInv_params.p_ref = 1.5 * (v_od_eq * i_od_eq + v_oq_eq * i_oq_eq);
GFMInv_params.q_ref = 1.5 * (v_oq_eq * i_od_eq - v_od_eq * i_oq_eq);
GFMInv_params.v_ref = sqrt(v_od_eq^2 + v_oq_eq^2);

% plot settings
color_unctrl = [0.0000, 0.0000, 0.6000];  % blue
color_ctrl   = [0.8500, 0.3250, 0.0980];  % orange
line_width   = 1.5;
font_size    = 12;
font_name    = 'Times New Roman';

%% Approximated transfer function model

% current loop
G_il = (GFMInv_params.kip*s + GFMInv_params.kii) / (GFMInv_params.Lf*s^2 + (GFMInv_params.kip+GFMInv_params.Rf)*s + GFMInv_params.kii);

% voltage loop
G_vl = G_il*(GFMInv_params.kvp*s + GFMInv_params.kvi) / (GFMInv_params.Cf*s^2 + G_il*GFMInv_params.kvp*s + G_il*GFMInv_params.kvi);

no_ctrl_tf = G_vl;

%% parameter tuning
% leaky-PI strategy
sigma_target = 5; 
freq_target = 100.0 * 2*pi;
param_a = 1;
param_b = 0.5;

[kp, ki, is_feasible] = func_solve_PI_passivity(no_ctrl_tf, sigma_target, freq_target, param_a, param_b);

if is_feasible
    func_plot_PI_feasible_region(no_ctrl_tf, sigma_target, freq_target, param_a, param_b);
end

%% Frequency-domain comparison

% construct PI controller
ctrler_tf = kp + ki / (param_a*s + param_b);

% open-loop transfer function with PI controller
ctrl_tf = ctrler_tf * no_ctrl_tf;

% Nyquist diagram comparison
fig_nyq = figure('Name', 'Nyquist Comparison');

nyquist(no_ctrl_tf, w_nyq); hold on;
lines_nyq1 = findobj(gca, 'Type', 'line');
set(lines_nyq1, 'Color', color_unctrl, 'LineWidth', line_width);

nyquist(ctrl_tf, w_nyq);
lines_nyq2 = setdiff(findobj(gca, 'Type', 'line'), lines_nyq1);
set(lines_nyq2, 'Color', color_ctrl, 'LineWidth', line_width);

grid off;
title('Nyquist Diagram Comparison', 'FontSize', font_size, 'FontName', font_name);
set(gca, 'FontSize', font_size, 'FontName', font_name);

% Positive Damping Region
theta_pd = linspace(0, 2*pi, 100); % fewer points for the circle
r_pd = 1 / (2 * sigma_target);
x_pd = r_pd + r_pd * cos(theta_pd);
y_pd = r_pd * sin(theta_pd);
h_pd = fill(x_pd, y_pd, [0.85, 0.90, 1.00], 'EdgeColor', 'none');
uistack(h_pd, 'bottom');

lgd = legend('w/o ctrl', 'with ctrl', 'PD Region', 'Location', 'best');
set(lgd, 'FontSize', font_size, 'FontName', font_name, 'Interpreter', 'tex');
xlim([-0.1, 1.2])

% Bode diagram comparison
figure('Name', 'Bode Comparison');

bode(no_ctrl_tf); hold on;
lines_bode1 = findobj(gcf, 'Type', 'line');
set(lines_bode1, 'Color', color_unctrl, 'LineWidth', line_width);

bode(ctrl_tf);
lines_bode2 = setdiff(findobj(gcf, 'Type', 'line'), lines_bode1);
set(lines_bode2, 'Color', color_ctrl, 'LineWidth', line_width);

lgd2 = legend('w/o ctrl', 'with ctrl', 'Location', 'best');
set(lgd2, 'FontSize', font_size, 'FontName', font_name, 'Interpreter', 'tex');
grid off;
title('Bode Diagram Comparison', 'FontSize', font_size, 'FontName', font_name);
set(findall(gcf,'type','axes'), 'FontSize', font_size, 'FontName', font_name);

%% Transient simulation
power_filter_tf = 1;
G_VQ_seq = [0, power_filter_tf * ctrler_tf];

legend_names = {'w/o ctrl', 'with ctrl'};
color_seq = {color_unctrl, color_ctrl};

figure('Name', 'Transient Response'); hold on;
ylabel('Angle (deg)', 'FontSize', font_size, 'FontName', font_name); grid off;
set(gca, 'FontSize', font_size, 'FontName', font_name);

for rr = 1:length(G_VQ_seq)
    % controller initialization
    if isnumeric(G_VQ_seq(rr))
        G_VQ_seq(rr) = tf(G_VQ_seq(rr));
    end
    [num, den] = tfdata(G_VQ_seq(rr), 'v');
    [A_d, B_d, C_d, D_d] = tf2ss(num, den);
    GFMInv_params.droop_A = A_d;
    GFMInv_params.droop_B = B_d;
    GFMInv_params.droop_C = C_d;
    GFMInv_params.droop_D = D_d;
    
    % re-compute equilibrium since dimensions change with different controllers
    x_equi = calc_equilibrium(GFMInv_params, System_params);
    
    % simulation
    [t, x] = ode15s(@(t, x) calc_dynamics(t, x, GFMInv_params, System_params), tspan, x_equi, opts);
    V_DQ = sqrt(x(:,1).^2 + x(:,2).^2);    % DQ voltage magnitude
    
    % plot
    plot(t, V_DQ, 'LineWidth', line_width, 'Color', color_seq{rr}, 'DisplayName', legend_names{rr});
    
    % Save uncontrolled system data for FFT analysis later
    if rr == 1
        t_unctrl = t;
        V_DQ_unctrl = V_DQ;
    end
end

lgd3 = legend('Location', 'best'); box on;
set(lgd3, 'FontSize', font_size, 'FontName', font_name);

%% Fourier analysis of transient response

% Resample the non-uniform ode15s time-steps to a uniform time step based on dt
Fs = 1/dt; 
t_uniform = tspan(1):dt:tspan(end);
V_DQ_uniform = interp1(t_unctrl, V_DQ_unctrl, t_uniform, 'linear');

% Remove DC component to clearly observe the oscillation frequency peaks
V_DQ_ac = V_DQ_uniform - mean(V_DQ_uniform);

% Perform Fast Fourier Transform (FFT)
L_fft = length(V_DQ_ac);
Y_fft = fft(V_DQ_ac);

% Compute the two-sided spectrum P2 and then single-sided spectrum P1
P2 = abs(Y_fft / L_fft);
P1 = P2(1:floor(L_fft/2)+1);
P1(2:end-1) = 2*P1(2:end-1);

% Define frequency axis
f_axis = Fs * (0:(floor(L_fft/2))) / L_fft;

% Plot Spectrum
figure('Name', 'FFT Analysis');
plot(f_axis, P1, 'Color', color_unctrl, 'LineWidth', line_width);
title('Frequency Spectrum of Uncontrolled Voltage', 'FontSize', font_size, 'FontName', font_name);
xlabel('Frequency (Hz)', 'FontSize', font_size, 'FontName', font_name);
ylabel('|V(f)| Amplitude', 'FontSize', font_size, 'FontName', font_name);
xlim([0, 200]); 
grid off; box on;
set(gca, 'FontSize', font_size, 'FontName', font_name);
