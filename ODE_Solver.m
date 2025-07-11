% High-order ODE solver via cascaded MRRs

close all; clear; clc;

%% ---------------- User parameters ----------------
order    = 1;                  % 1, 2, or 3
k        = [0.5, 0.3, 0.2];    % [ns^-1] per stage
A        = 1e9;                % ns → s scaling
r_tolerances = [0]; % ± coupling tolerance (%)
loss_tolerance =0.01;               % ± intrinsic-loss tolerance (%)
detuning_tolerance  = 0;                % detuning (%)
N_monte_carlo = 1000;               % number of Monte Carlo trials
R     = [5e-3, 4e-3, 3e-3];    % ring radii [m]
neff  = 1.5;                   % effective index

%% ---------------- Input signal -------------------
C        = 40;                  % step amplitude
%x_fun = @(t) C * (t > 0);      % unit-step of amplitude C at t=0
f0 = 2.5e9;                     % sine frequency [Hz]
x_fun = @(t) C * sin(2*pi*f0*t);

%% ---------------- Static constants ---------
c     = 3e8;                   % speed of light [m/s]

%% ---------------- MRR geometry & derived ---------
L     = 2*pi*R;                         % round-trip length [m]
k_i   = k * A;                          % [s^-1]
tau_c = 1 ./ k_i;                       % cavity lifetime [s]
tau_rt= L ./ (c/neff);                  % round-trip time [s]
tau_n = tau_c ./ tau_rt;                % normalized lifetime
r     = sqrt(tau_n ./ (1 + tau_n));     % nominal coupling coeff
alpha = ones(size(r));                  % nominal intrinsic loss

%% → compute ordinal suffix for pretty printing
suffix = {'st','nd','rd','th'};
if mod(order,100)>=11 && mod(order,100)<=13
    suff = 'th';
else
    suff = suffix{min(mod(order,10)+1,4)};
end

%% ---------------- Simulation grid ----------------
N    = 1e5;
tmin = -100e-9;  tmax = 100e-9;
time = linspace(tmin, tmax, N);
dt   = time(2) - time(1);
Df   = linspace(-1/(2*dt), 1/(2*dt), N);

% Evaluate input signal
in_t = x_fun(time);
% FFT of input signal
X   = fftshift(fft(in_t));

P = abs(X).^2;        % Power spectrum

% Normalize to peak power
P = P / max(P);

% Convert to dB scale
P_dB = 10 * log10(P);

% Find frequencies above -3 dB threshold
threshold_dB = -3;
mask = P_dB >= threshold_dB;

% Extract frequency range above -3 dB
f_bw = Df(mask);
bandwidth_Hz = max(f_bw) - min(f_bw);

%% ---------------- State-space ODE45 --------------
ss_fun = @(t,y) ss_odes(t, y, k(1:order), A, x_fun(t));
y0     = zeros(order,1);
[t_ode, y_state] = ode45(ss_fun, [tmin tmax], y0);
y_ode45 = interp1(t_ode, y_state(:,end), time);

%% ---------------- Ideal ODE TF -------------------
H_ode = ones(1,N);
for i = 1:order
    H_ode = H_ode .* ((1/k(i)) .* (1/tau_c(i) ./ (1/tau_c(i) + 1j*2*pi*Df)));
end
y_ode_tf = real(ifft(fftshift(X .* H_ode)));

%% ---------------- Cascaded MRR TF ----------------
beta  = 2*pi*Df / (c/neff);
H_mrr = ones(1,N);
for i = 1:order
    H_drop = (1/k(i)) .* ((1-r(i)^2).*alpha(i) ./ ...
             (1 - r(i)^2.*alpha(i).*exp(-1j*beta*L(i))));
    H_mrr  = H_mrr .* H_drop;
end
y_mrr = real(ifft(fftshift(X .* H_mrr)));

%% ---------------- Nominal plots ------------------
% Normalize and convert to dB
[H_ode_Df_GHz, H_ode_dB] = normalize_power_db(Df, H_ode);
[H_mrr_Df_GHz, H_mrr_dB] = normalize_power_db(Df, H_mrr);

figure('Name', 'Frequency-Domain');

% Plot ideal ODE-based transfer function
plot(H_ode_Df_GHz, H_ode_dB, 'b', 'LineWidth', 1.4); 
hold on;

% Plot cascaded MRR transfer function
plot(H_mrr_Df_GHz, H_mrr_dB, 'r', 'LineWidth', 1.4);

%Plot the signal bandwidth
plot(Df / 1e9, P_dB, 'g', 'LineWidth', 1.4);
yline(-3, 'r--', '-3 dB', 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom');
% Bandwidth limits
min_f = min(f_bw);
max_f = max(f_bw);
% Compute bandwidth in GHz
bw_GHz = (max_f - min_f)/1e9;

% Compute center of the shaded area
center_f = (min_f + max_f)/2 / 1e9;  % in GHz
text_y = -120;                         % vertical position for the label

% Add text label
text(center_f, text_y, ['B = ' num2str(bw_GHz, '%.2f') ' GHz'], ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
fill_x = [min_f max_f max_f min_f]/1e9;       % GHz for plotting
fill_y = [-200 -200 0 0];
fill(fill_x, fill_y, [0.2 0.4 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
xline(min(f_bw)/1e9, 'k--');
xline(max(f_bw)/1e9, 'k--');

% Axis labels and styling
xlabel('Frequency [GHz]');
ylabel('Magnitude [dB]');

title('Frequency-Domain Transfer Functions');
legend({ ...
    'Ideal ODE TF', ...
    'Cascaded MRR TF', ...
    'Input Signal Spectrum', ...
    '-3 dB Reference', ...
    '3 dB Bandwidth Region' ...
}, 'Location', 'best');
grid on;

xlim([-15, 15]);
ylim([-30, 0]);

figure('Name','Time-Domain Signals');

subplot(3,1,1);
  plot(time*1e9, in_t,'k','LineWidth',2);
  ylabel('x(t)'); xlim([0 tmax*1e9]); grid on;
subplot(3,1,2);
  plot(time*1e9, y_ode45,'k','LineWidth',1.5); hold on;
  plot(time*1e9, y_mrr,  'r','LineWidth',1.5);
  ylabel('y(t)'); legend('ODE45','MRR'); xlim([0 tmax*1e9]); grid on;
subplot(3,1,3);
  plot(time*1e9, abs(y_ode45).^2,'k','LineWidth',1.5); hold on;
  plot(time*1e9, abs(y_mrr).^2,  'r','LineWidth',1.5);
  ylabel('|y|^2'); xlabel('Time [ns]');
  legend('ODE45','MRR'); xlim([0 tmax*1e9]); grid on;

%% ---------------- Monte-Carlo analysis -----------
gain_ideal = C / prod(k(1:order));
% Preallocate results
gain_monte_carlo_all = cell(length(r_tolerances),1);
H_all_all      = cell(length(r_tolerances),1);
y_all_all      = cell(length(r_tolerances),1);

for idx = 1:length(r_tolerances)
    tol_r = r_tolerances(idx);
    [gain_monte_carlo, H_all, y_all, H_env_max, H_env_min, y_env_max, y_env_min] = ...
        MonteCarlo(N_monte_carlo, k, order, time, N, r, tol_r, alpha, loss_tolerance, neff, detuning_tolerance, Df, c, L, X, H_mrr);
    gain_monte_carlo_all{idx} = gain_monte_carlo;
    H_all_all{idx}      = H_all;
    y_all_all{idx}      = y_all;
    fprintf('Monte-Carlo run %d/%d with tol_r = %.3f done.\n', idx, length(r_tolerances), tol_r);

    % Error statistics
    mean_err = mean((gain_monte_carlo - gain_ideal)/gain_ideal)*100;
    std_err  = std((gain_monte_carlo - gain_ideal)/gain_ideal)*100;
    fprintf('\nMonte-Carlo (%d trials, %d^%s-order):\n', N_monte_carlo, order, suff);
    fprintf('   Ideal DC gain = %.3e\n', gain_ideal);
    fprintf('   Mean error    = %+6.2f %%\n', mean_err);
    fprintf('   Std  error    =  ±%.2f %%\n\n', std_err);
end

for idx = 1:length(r_tolerances)
    %% ---------------- Percentile envelopes ----------------
    p_low  =  5;   p_mid = 50;   p_high = 95;
    YL = prctile(y_all_all{idx}, p_low,  1);
    YM = prctile(y_all_all{idx}, p_mid,  1);
    YH = prctile(y_all_all{idx}, p_high, 1);
    FL = prctile(H_all_all{idx}, p_low,  1);
    FM = prctile(H_all_all{idx}, p_mid,  1);
    FH = prctile(H_all_all{idx}, p_high, 1);
    
    %% ---------------- Plot percentile envelopes ------------
    figure('Name', sprintf('Percentile Envelopes & Density - Coupling Tolerance: %.2f%%', r_tolerances(idx)*100));
    
    % Frequency-domain envelope
    subplot(3,1,1);
    [H_ode_Df_GHz, H_ode_dB] = normalize_power_db(Df, H_ode);
    [H_mrr_Df_GHz, H_mrr_dB] = normalize_power_db(Df, H_mrr);
    plot(H_ode_Df_GHz, H_ode_dB,'b','LineWidth',1.4); hold on;
    plot(H_mrr_Df_GHz, H_mrr_dB, 'r','LineWidth',1.4); hold on;
    
    fill([Df/1e9, fliplr(Df/1e9)], [10*log10(FH), fliplr(10*log10(FL))], ...
         'r','FaceAlpha',0.2,'EdgeColor','none');
    plot(Df/1e9, 10*log10(FM), 'r--', 'LineWidth',1);
    xlabel('Frequency [GHz]'); ylabel('Magnitude [dB]');
    legend('Ideal ODE TF','Nominal','5–95% envelope','Median','Location','SouthWest');
    xlim([-15 15]); ylim([-30 0]); grid on;
    
    % Time-domain envelope
    subplot(3,1,2);
    plot(time*1e9, y_ode45, 'k','LineWidth',1.5); hold on;
    plot(time*1e9, YM,      'r--','LineWidth',1);
    fill([time*1e9, fliplr(time*1e9)], [YH, fliplr(YL)], ...
         'r','FaceAlpha',0.2,'EdgeColor','none');
    xlabel('Time [ns]'); ylabel('y(t)');
    legend('ODE45','Median MC','5–95% envelope','Location','SouthEast');
    xlim([0 tmax*1e9]); grid on;
    
    % 2D density heat-map (time vs. amplitude)
    subplot(3,1,3);
    edges = linspace(min(y_all_all{idx}(:)),max(y_all_all{idx}(:)),100);
    pd = zeros(length(edges)-1, N);
    for i = 1:N
        counts = histcounts(y_all_all{idx}(:,i), edges);
        pd(:,i) = counts / N;
    end
    imagesc(time*1e9, edges(1:end-1), pd);
    axis xy; hold on;
    plot(time*1e9, y_ode45, 'w--','LineWidth',1.5);  % ODE45 guideline
    xlabel('Time [ns]'); ylabel('y amplitude');
    title('Probability density of y(t) with ODE45 guideline');
    xlim([0 tmax*1e9]); grid on;
    colorbar;
end

%% ---------------- Plot different montecarlo analysis ----------------
figure('Name','Mean error vs. coupling tolerance');

mse_errors = zeros(length(r_tolerances), 1);
rmsd = zeros(length(r_tolerances), 1);
ci_half_widths = zeros(length(r_tolerances), 1);

t_alpha = 0.05;  % 95% confidence level

for idx = 1:length(r_tolerances)
    % ⚠️ Compute relative error in percentage
    rel_err_percent = (gain_monte_carlo_all{idx} - gain_ideal) / gain_ideal * 100;

    % Mean Squared Error (in %²)
    mse_errors(idx) = mean(rel_err_percent.^2);

    % Root Mean Squared Error (in %)
    rmsd(idx) = sqrt(mse_errors(idx));

    % Standard error of the squared relative error
    N = length(rel_err_percent);
    sem = std(sqrt(rel_err_percent.^2)) / sqrt(N);

    % t-value for confidence interval
    tval = tinv(1 - t_alpha/2, N - 1);

    % Half-width of CI for MSE
    ci_half_widths(idx) = tval * sem;
end

% Plot with error bars
errorbar(r_tolerances *100, rmsd, ci_half_widths, 'o-', ...
    'LineWidth', 1.5, 'MarkerSize', 6, 'Color', [0.2 0.6 0.8]);

xlabel('Coupling Tolerance (%)');
ylabel('Root Mean Square Error (%)');
title('Root Mean Square Error vs. Coupling Tolerance');
grid on;

%% ---------------- helper function -----------------------------------
function dydt = ss_odes(~, y, kvec, A, xin)
    N    = numel(kvec);
    dydt = zeros(N,1);
    dydt(1) = A*( xin       - kvec(1)*y(1) );
    for j = 2:N
        dydt(j) = A*( y(j-1)   - kvec(j)*y(j) );
    end
end

function [f_GHz, power_dB] = normalize_power_db(freq_Hz, H)
%NORMALIZE_POWER_DB Normalizes a frequency response and converts it to dB
%
%   Inputs:
%       freq_Hz : frequency vector in Hz
%       H       : complex frequency-domain transfer function (same size)
%
%   Outputs:
%       f_GHz    : frequency vector in GHz
%       power_dB : normalized power response in dB

    % Convert frequency to GHz
    f_GHz = freq_Hz / 1e9;

    % Normalize magnitude response
    H_norm = abs(H) / max(abs(H));

    % Convert to power in dB
    power_dB = 10 * log10(H_norm .^ 2);
end