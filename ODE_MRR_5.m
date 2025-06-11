%% ODE_MRR_highorder_fixed2.m
% Solve up to 3rd-order ODE via ODE45 (state‐space) vs. cascaded MRR TF
% Uses the same 1/k prefactors in both.

close all; clear; clc;

%% ---------------- User parameters ----------------
order    = 3;                % 1, 2 or 3
k        = [0.5, 0.3, 0.2];  % [ns^-1] per stage
A        = 1e9;              % ns → s
C        = 4;                % Step amplitude
tol_r    = 0.05;             % ± coupling tolerance
tol_loss = 0.05;             % ± loss tolerance
tol_det  = 1e8;              % ± detuning [Hz]

%% ---------------- Input signal ----------------
x_fun = @(t) C*(t>0);        % Step at t=0

%% ---------------- MRR geometry & derived ----------------
R         = [5e-3,4e-3,3e-3];    % m
neff      = 1.5;
c         = 3e8;
L         = 2*pi*R;              % round‐trip length [m]
k_i       = k * A;               % [s^-1]
tau_c     = 1./(k_i);            % cavity lifetime [s]
tau_rt    = L./(c/neff);         % round‐trip time [s]
tau_n     = tau_c ./ tau_rt;
r_nom     = sqrt(tau_n./(1+tau_n));
alpha_nom = ones(size(r_nom));
detune_nom= zeros(size(r_nom));

%% ---------------- Simulation grid ----------------
N    = 1e5;
tmin = -100e-9; tmax = 100e-9;
time = linspace(tmin, tmax, N);
dt   = time(2)-time(1);
Df   = linspace(-1/(2*dt), 1/(2*dt), N);

% Input in time & frequency
in_t = x_fun(time);
IN   = fftshift(fft(in_t));

%% ---------------- State‐space ODE (ODE45) ----------------
% Build cascaded ODEs EXACTLY like professor's: A*(input - k*y)
ss_fun = @(t,y) ss_odes(t, y, k(1:order), A, x_fun(t));

% initial state all zero
y0 = zeros(order,1);

% solve on [tmin tmax] in seconds
[t_ode, y_state] = ode45(ss_fun, [tmin tmax], y0);

% final output is the last state
y_ode45 = interp1(t_ode, y_state(:,end), time);

%% ---------------- Ideal ODE TF (for reference) ----------------
H_ode = ones(1,N);
for i=1:order
  H_ode = H_ode .* (1/k(i))*(1/tau_c(i))./(1/tau_c(i)+1j*2*pi*Df);
end
y_ode_tf = real(ifft(fftshift(IN .* H_ode)));

%% ---------------- Cascaded MRR TF ----------------
beta  = 2*pi*Df/(c/neff);
H_mrr = ones(1,N);
for i=1:order
  H_drop = (1/k(i)) * ...
           ((1-r_nom(i)^2)*alpha_nom(i) ./ ...
            (1 - r_nom(i)^2*alpha_nom(i).*exp(-1j*beta*L(i))));
  H_mrr = H_mrr .* H_drop;
end
y_mrr = real(ifft(fftshift(IN .* H_mrr)));

%% ---------------- Plots ----------------

% 1) Nominal comparison ODE45 vs MRR
figure('Name','Nominal Response');
plot(time*1e9, y_ode45,    'k','LineWidth',2); hold on;
plot(time*1e9, y_mrr,      'r--','LineWidth',2);
xlabel('Time [ns]'); ylabel('y(t)');
legend('ODE45','MRR','Location','Southeast');

% Determine ordinal suffix
if mod(order, 100) >= 11 && mod(order, 100) <= 13
    suffix = 'th';
else
    switch mod(order, 10)
        case 1
            suffix = 'st';
        case 2
            suffix = 'nd';
        case 3
            suffix = 'rd';
        otherwise
            suffix = 'th';
    end
end
title(sprintf('%d^{%s}-order: ODE45 vs. MRR', order, suffix));
grid on;

% 2) Ideal‐TF vs MRR TF
figure('Name','Frequency-Domain TF');
plot(Df/1e9, 10*log10(abs(H_ode./max(abs(H_ode))).^2),'b','LineWidth',1.5); hold on;
plot(Df/1e9, 10*log10(abs(H_mrr./max(abs(H_mrr))).^2),'r','LineWidth',1.5);
xlabel('Freq [GHz]'); ylabel('Magnitude [dB]');
legend('Ideal ODE TF','Cascaded MRR TF','Location','Southeast');
grid on; xlim([-15 15]); ylim([-30 0]);

% 3) Time‐domain input + squared output
figure('Name','Time-Domain Signals');
subplot(3,1,1);
  plot(time*1e9,in_t,'k','LineWidth',2);
  ylabel('x(t)'); xlim([0 tmax*1e9]); grid on;
subplot(3,1,2);
  plot(time*1e9,y_ode45,'k','LineWidth',1.5); hold on;
  plot(time*1e9,y_mrr,  'r--','LineWidth',1.5);
  ylabel('y(t)'); legend('ODE45','MRR'); xlim([0 tmax*1e9]); grid on;
subplot(3,1,3);
  plot(time*1e9,abs(y_ode45).^2,'k','LineWidth',1.5); hold on;
  plot(time*1e9,abs(y_mrr).^2,    'r--','LineWidth',1.5);
  ylabel('|y|^2'); xlabel('Time [ns]'); legend('ODE45','MRR'); xlim([0 tmax*1e9]); grid on;

%% ---------------- (Optional) Monte Carlo ----------------
% — as before, but now both TFs share same 1/k prefactors —  

%% ---------------- State‐space helper ----------------
function dydt = ss_odes(~, y, kvec, A, xin)
  N    = numel(kvec);
  dydt = zeros(N,1);
  % first stage
  dydt(1) = A*( xin - kvec(1)*y(1) );
  % subsequent stages
  for j=2:N
    dydt(j) = A*( y(j-1) - kvec(j)*y(j) );
  end
end
