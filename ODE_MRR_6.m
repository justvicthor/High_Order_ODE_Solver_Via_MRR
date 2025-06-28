%% ODE_MRR_highorder_full_CorrectedDC.m
% High-order ODE solver via cascaded MRRs
% Steady-state DC gain now includes C
% + Monte-Carlo  + Parametric sweeps + Active tuning + Speed benchmark

close all; clear; clc;

%% ---------------- User parameters ----------------
order    = 3;                  % 1, 2, or 3
k        = [0.5, 0.3, 0.2];    % [ns^-1]  per stage
A        = 1e9;                % ns → s scaling
C        = 4;                  % step amplitude
tol_r    = 0.05;               % ± coupling tolerance (relative)
tol_loss = 0.05;               % ± intrinsic-loss tolerance (relative)
tol_det  = 1e8;                % detuning σ [Hz]

%% ---------------- Input signal -------------------
x_fun = @(t) C * (t > 0);      % unit-step of amplitude C at t=0

%% ---------------- MRR geometry & derived ---------
R     = [5e-3, 4e-3, 3e-3];    % ring radii [m]
neff  = 1.5;                   % effective index
c     = 3e8;                   % speed of light [m/s]
L     = 2*pi*R;                % round-trip length [m]

k_i   = k * A;                 % [s^-1]
tau_c = 1 ./ k_i;              % cavity lifetime [s]
tau_rt= L ./ (c/neff);         % round-trip time [s]
tau_n = tau_c ./ tau_rt;       % normalized lifetime

r_nom     = sqrt(tau_n ./ (1 + tau_n));   % nominal coupling coeff
alpha_nom = ones(size(r_nom));            % nominal intrinsic loss
detune_nom= zeros(size(r_nom));           % nominal detuning

%% ---------------- Simulation grid ----------------
N    = 1e5;
tmin = -100e-9;  tmax = 100e-9;
time = linspace(tmin, tmax, N);
dt   = time(2) - time(1);
Df   = linspace(-1/(2*dt), 1/(2*dt), N);

in_t = x_fun(time);
IN   = fftshift(fft(in_t));               % FFT of input

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
y_ode_tf = real(ifft(fftshift(IN .* H_ode)));

%% ---------------- Cascaded MRR TF ----------------
beta  = 2*pi*Df / (c/neff);
H_mrr = ones(1,N);
for i = 1:order
    H_drop = (1/k(i)) .* ((1-r_nom(i)^2).*alpha_nom(i) ./ ...
             (1 - r_nom(i)^2.*alpha_nom(i).*exp(-1j*beta*L(i))));
    H_mrr  = H_mrr .* H_drop;
end
y_mrr = real(ifft(fftshift(IN .* H_mrr)));

%% ---------------- Nominal plots ------------------
figure('Name','Nominal Response');
plot(time*1e9, y_ode45,'k','LineWidth',2); hold on;
plot(time*1e9, y_mrr,  'r--','LineWidth',2);
xlabel('Time [ns]'); ylabel('y(t)');
legend('ODE45','MRR','Location','Southeast');
suffix = {'st','nd','rd','th'}; 
if mod(order,100)>=11 && mod(order,100)<=13, suff='th'; 
else suff=suffix{min(mod(order,10),4)}; end
title(sprintf('%d^{%s}-order: ODE45 vs MRR',order,suff)); grid on;

figure('Name','Frequency-Domain TF');
plot(Df/1e9,10*log10(abs(H_ode/max(abs(H_ode))).^2),'b','LineWidth',1.4); hold on;
plot(Df/1e9,10*log10(abs(H_mrr/max(abs(H_mrr))).^2),'r','LineWidth',1.4);
xlabel('Freq [GHz]'); ylabel('Magnitude [dB]');
legend('Ideal ODE TF','Cascaded MRR TF'); grid on; xlim([-15 15]); ylim([-30 0]);

%% ---------------- Time-domain Signals ----------------
figure('Name','Time-Domain Signals');
subplot(3,1,1);
  plot(time*1e9, in_t,'k','LineWidth',2);
  ylabel('x(t)'); xlim([0 tmax*1e9]); grid on;
subplot(3,1,2);
  plot(time*1e9, y_ode45,'k','LineWidth',1.5); hold on;
  plot(time*1e9, y_mrr,  'r--','LineWidth',1.5);
  ylabel('y(t)'); legend('ODE45','MRR'); xlim([0 tmax*1e9]); grid on;
subplot(3,1,3);
  plot(time*1e9, abs(y_ode45).^2,'k','LineWidth',1.5); hold on;
  plot(time*1e9, abs(y_mrr).^2,  'r--','LineWidth',1.5);
  ylabel('|y|^2'); xlabel('Time [ns]');
  legend('ODE45','MRR'); xlim([0 tmax*1e9]); grid on;

%% ---------------- Monte-Carlo analysis -----------
Nmc = 100000;
% Ideal DC gain now includes C:
dc_gain_ideal = C / prod(k(1:order));
dc_gain_mc    = zeros(Nmc,1);
beta = 2*pi*Df/(c/neff);
idx_dc = find(time>=1e-9 & time<=5e-9);

for mc = 1:Nmc
    r_mc      = r_nom(1:order)   .* (1 + tol_r   * randn(1,order));
    alpha_mc  = alpha_nom(1:order).* (1 + tol_loss* randn(1,order));
    detune_mc = tol_det * randn(1,order);
    
    H_mc = ones(1,N);
    for i = 1:order
        beta_i = 2*pi*(Df + detune_mc(i))/(c/neff);
        H_comp = (1/k(i)) .* ((1-r_mc(i)^2).*alpha_mc(i) ./ ...
                 (1 - r_mc(i)^2.*alpha_mc(i).*exp(-1j*beta_i*L(i))));
        H_mc   = H_mc .* H_comp;
    end
    y_mc = real(ifft(fftshift(IN .* H_mc)));
    % No division by C here—y_mc is in absolute units
    dc_gain_mc(mc) = mean(y_mc(idx_dc));
end

mean_err = mean((dc_gain_mc - dc_gain_ideal)/dc_gain_ideal)*100;
std_err  = std((dc_gain_mc - dc_gain_ideal)/dc_gain_ideal)*100;
fprintf('\nMonte-Carlo (%d trials, %d^%s-order):\n',Nmc,order,suff);
fprintf('   Ideal DC gain = %.3e\n',dc_gain_ideal);
fprintf('   Mean error    = %+6.2f %%\n',mean_err);
fprintf('   Std  error    =  ±%.2f %%\n\n',std_err);

figure('Name','DC-gain error histogram');
histogram((dc_gain_mc - dc_gain_ideal)/dc_gain_ideal*100,50);
xlabel('DC-gain error (%)'); ylabel('# trials'); grid on;

%% ---------------------------------------------------------------------
%% EXTRA STUDIES  – toggle below
do_parametric_sweep = false;
do_active_tuning    = false;
do_speed_benchmark  = false;
%% ---------------------------------------------------------------------

%% 1) PARAMETRIC SWEEPS -------------------------------------------------
if do_parametric_sweep
    fprintf('=== PARAMETRIC SWEEP (order = %d) ===\n',order);
    
    sweep_r    = linspace(-0.10,+0.10,11);
    sweep_loss = linspace(-0.10,+0.10,11);
    sweep_det  = linspace(-3e8,+3e8,13);
    
    err_r = zeros(size(sweep_r));
    err_l = zeros(size(sweep_loss));
    err_d = zeros(size(sweep_det));
    
    for ii = 1:numel(sweep_r)
        r_pert = r_nom(1:order).*(1 + sweep_r(ii));
        Htmp   = ones(1,N);
        for i=1:order
            H_comp = (1/k(i)) .* ((1-r_pert(i)^2) ./ ...
                     (1 - r_pert(i)^2.*exp(-1j*beta*L(i))));
            Htmp   = Htmp .* H_comp;
        end
        ytmp = real(ifft(fftshift(IN .* Htmp)));
        err_r(ii) = (mean(ytmp(idx_dc)) - dc_gain_ideal)/dc_gain_ideal*100;
    end
    
    for ii = 1:numel(sweep_loss)
        a_pert = alpha_nom(1:order).*(1 + sweep_loss(ii));
        Htmp   = ones(1,N);
        for i=1:order
            H_comp = (1/k(i)) .* ((1-r_nom(i)^2).*a_pert(i) ./ ...
                     (1 - r_nom(i)^2.*a_pert(i).*exp(-1j*beta*L(i))));
            Htmp   = Htmp .* H_comp;
        end
        ytmp = real(ifft(fftshift(IN .* Htmp)));
        err_l(ii) = (mean(ytmp(idx_dc)) - dc_gain_ideal)/dc_gain_ideal*100;
    end
    
    for ii = 1:numel(sweep_det)
        Htmp = ones(1,N);
        for i=1:order
            beta_i = 2*pi*(Df + sweep_det(ii))/(c/neff);
            H_comp = (1/k(i)) .* ((1-r_nom(i)^2) ./ ...
                     (1 - r_nom(i)^2.*exp(-1j*beta_i*L(i))));
            Htmp   = Htmp .* H_comp;
        end
        ytmp = real(ifft(fftshift(IN .* Htmp)));
        err_d(ii) = (mean(ytmp(idx_dc)) - dc_gain_ideal)/dc_gain_ideal*100;
    end
    
    figure('Name','Parametric sweeps');
    subplot(1,3,1); plot(sweep_r*100,err_r,'o-'); grid on;
      xlabel('Δr (%)'); ylabel('DC-gain error (%)'); title('Coupling');
    subplot(1,3,2); plot(sweep_loss*100,err_l,'o-'); grid on;
      xlabel('Δα (%)'); title('Loss');
    subplot(1,3,3); plot(sweep_det/1e6,err_d,'o-'); grid on;
      xlabel('Δf (MHz)'); title('Detuning');
end

%% 2) ACTIVE-TUNING STUDY ----------------------------------------------
if do_active_tuning
    fprintf('=== ACTIVE-TUNING STUDY (order = %d) ===\n',order);
    delta_f_max = 1.5e8;   % ±150 MHz tuning
    Nmc_tune    = 500;
    dc_raw   = zeros(Nmc_tune,1);
    dc_tuned = zeros(Nmc_tune,1);

    for mc = 1:Nmc_tune
        r_mc     = r_nom(1:order).*(1 + tol_r*randn(1,order));
        alpha_mc = alpha_nom(1:order).*(1 + tol_loss*randn(1,order));
        det_mc   = tol_det*randn(1,order);

        % raw
        Hraw = ones(1,N);
        for i=1:order
            beta_i = 2*pi*(Df + det_mc(i))/(c/neff);
            H_comp = (1/k(i)) .* ((1-r_mc(i)^2).*alpha_mc(i) ./ ...
                     (1 - r_mc(i)^2.*alpha_mc(i).*exp(-1j*beta_i*L(i))));
            Hraw   = Hraw .* H_comp;
        end
        yraw = real(ifft(fftshift(IN .* Hraw)));
        dc_raw(mc) = mean(yraw(idx_dc));

        % tuned
        tune_cmd = min(max(-det_mc,-delta_f_max),+delta_f_max);
        Htun = ones(1,N);
        for i=1:order
            beta_i = 2*pi*(Df + det_mc(i)+tune_cmd(i))/(c/neff);
            H_comp = (1/k(i)) .* ((1-r_mc(i)^2).*alpha_mc(i) ./ ...
                     (1 - r_mc(i)^2.*alpha_mc(i).*exp(-1j*beta_i*L(i))));
            Htun   = Htun .* H_comp;
        end
        ytun = real(ifft(fftshift(IN .* Htun)));
        dc_tuned(mc) = mean(ytun(idx_dc));
    end

    fprintf('   Mean raw   = %.2f  (err %+5.1f %%)\n', mean(dc_raw),   (mean(dc_raw)-dc_gain_ideal)/dc_gain_ideal*100);
    fprintf('   Mean tuned = %.2f  (err %+5.1f %%)\n\n', mean(dc_tuned), (mean(dc_tuned)-dc_gain_ideal)/dc_gain_ideal*100);

    figure('Name','Active tuning effect');
    histogram((dc_raw   - dc_gain_ideal)/dc_gain_ideal*100,40,'FaceColor','r','FaceAlpha',0.4); hold on;
    histogram((dc_tuned - dc_gain_ideal)/dc_gain_ideal*100,40,'FaceColor','g','FaceAlpha',0.4);
    legend('No tuning','With tuning'); xlabel('DC-gain error (%)'); grid on;
    title(sprintf('Tuning ±%.0f MHz',delta_f_max/1e6));
end

%% 3) SPEED BENCHMARK ---------------------------------------------------
if do_speed_benchmark
    fprintf('=== SPEED BENCHMARK ===\n');
    Nlong  = 2^19;
    time_l = linspace(tmin,tmax,Nlong);
    in_l   = x_fun(time_l);           IN_l = fftshift(fft(in_l));
    Df_l   = linspace(-1/(2*(time_l(2)-time_l(1))), 1/(2*(time_l(2)-time_l(1))), Nlong);
    beta_l = 2*pi*Df_l/(c/neff);

    tic;
    H_l = ones(1,Nlong);
    for i=1:order
        H_comp = (1/k(i)) .* ((1-r_nom(i)^2) ./ ...
                 (1 - r_nom(i)^2.*exp(-1j*beta_l*L(i))));
        H_l = H_l .* H_comp;
    end
    y_fft = real(ifft(fftshift(IN_l .* H_l)));
    t_fft = toc;

    tic;
    ss_fun_long = @(t,y) ss_odes(t,y,k(1:order),A,x_fun(t));
    [t_ode_l, y_state_l] = ode45(ss_fun_long,[tmin tmax],y0);
    y_ode_l = interp1(t_ode_l,y_state_l(:,end), time_l);
    t_ode45 = toc;

    fprintf('   FFT path: %.3f s  (N = %d)\n', t_fft, Nlong);
    fprintf('   ODE45  : %.3f s  (N = %d)\n\n', t_ode45, Nlong);
end

%% ---------------- helper function -----------------------------------
function dydt = ss_odes(~, y, kvec, A, xin)
    N    = numel(kvec);
    dydt = zeros(N,1);
    dydt(1) = A*( xin - kvec(1)*y(1) );
    for j = 2:N
        dydt(j) = A*( y(j-1) - kvec(j)*y(j) );
    end
end

