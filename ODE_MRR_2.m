%% ODE_MRR_highorder_with_basic_plots.m
% Extended MATLAB script to solve and analyze higher-order ODEs using cascaded MRRs
% Adds sensitivity to tolerances, intrinsic loss, and resonance detuning
% Also includes basic time-domain and frequency-domain plots similar to the professor's original code

close all; clear; clc;

%% ---------------- User parameters ----------------
order = 1;                 % Order of the ODE (1, 2, or 3)
k = [0.5, 0.3, 0.2];       % Coefficients (ns^-1) for each first-order stage
A = 1e9;                   % Time-scaling parameter (ns^-1)
C = 4;                     % Step input amplitude

tol_r    = 0.05;           % ± tolerance in coupling coefficient r (fraction)
tol_loss = 0.05;           % ± tolerance in round-trip loss alpha_i (fraction)
tol_det  = 1e8;            % ± tolerance in detuning (Hz)

%% ---------------- Input signal ----------------
% Define step input as function handle
x_fun = @(t) C*(t > 0);

%% ---------------- Derived ODE and MRR parameters ----------------
% Convert decay coefficients to absolute rates [s^-1]
k_i = k * A;

% MRR geometry & waveguide
R      = [5e-3, 4e-3, 3e-3];      % Ring radii [m]
neff   = 1.5;                     % Effective index
c      = 3e8;                     % Speed of light [m/s]
L      = 2*pi*R;                  % Round-trip length [m]

tau_c  = 1./(k_i);                % Cavity lifetime [s]
tau_rt = L./(c/neff);             % Round-trip time [s]

% Nominal coupling and loss
tau_n     = tau_c ./ tau_rt;r_nom     = sqrt(tau_n ./ (1 + tau_n));  % Nominal coupling coeff.
alpha_nom = ones(size(r_nom));           % Nominal intrinsic loss
detune_nom = zeros(size(R));             % Nominal resonance detuning [Hz]

%% ---------------- Simulation grid ----------------
N     = 1e5;
tmin  = -100e-9; tmax = 100e-9;
time  = linspace(tmin, tmax, N);
dt    = time(2) - time(1);
Df    = linspace(-1/(2*dt), 1/(2*dt), N);

% Input spectrum
in_t = x_fun(time);
IN   = fftshift(fft(in_t));

%% ---------------- Monte Carlo sensitivity ----------------
Nsamp = 50;                     % Number of Monte Carlo samples
y_ring_MC = zeros(N, Nsamp);
y_ode_MC  = zeros(N, Nsamp);

for m = 1:Nsamp
    % Randomize parameters within tolerance
    r_i     = r_nom(1:order) .* (1 + (2*rand(1,order)-1)*tol_r);
    alpha_i = (1 - tol_loss) + 2*tol_loss*rand(1,order);
    det_i   = detune_nom(1:order) + (2*rand(1,order)-1)*tol_det;

    % Build cascaded MRR TF
    H_tot = ones(1, N);
    for i = 1:order
        beta   = 2*pi*(Df - det_i(i)) ./ (c/neff);
        H_drop = (1 - r_i(i)^2) .* alpha_i(i) ./ (1 - r_i(i)^2 .* alpha_i(i) .* exp(-1j*beta*L(i)));
        H_tot  = H_tot .* H_drop;
    end

    % Ideal ODE TF for comparison
    H_ode = ones(1, N);
    for i = 1:order
        H_ode = H_ode .* (1/k(i)) * (1/tau_c(i)) ./ (1/tau_c(i) + 1j*2*pi*Df);
    end

    % Compute time-domain outputs
    y_ring_MC(:,m) = real(ifft(fftshift(IN .* H_tot)));
    y_ode_MC(:,m)  = real(ifft(fftshift(IN .* H_ode)));
end

%% ---------------- Nominal solution & ODE45 comparison ----------------
% Solve state-space ODE with ODE45 for nominal parameters
define_ss = @(t, y) ss_odes(t, y, k_i(1:order), x_fun(t));
%[t_ode, y_state] = ode45(define_ss, time, zeros(order,1));
y0 = zeros(order,1);
y0(1) = 1;  % Set initial condition for y(0)
[t_ode, y_state] = ode45(define_ss, time, y0);

y_ode45   = y_state(:, end);        % Last state as system output
y_ring_nom = y_ring_MC(:,1);        % Nominal MRR response

%% ---------------- Plot results: high-order comparison ----------------
figure('Name','High-Order ODE Solver');
subplot(2,1,1); hold on; grid on;
plot(time*1e9, y_ode45, 'k', 'LineWidth', 2);
plot(time*1e9, y_ring_nom, 'r--', 'LineWidth', 2);
legend('ODE45', 'MRR'); xlabel('Time [ns]'); ylabel('y(t)');
title(sprintf('%d^{th}-order ODE: Nominal Response', order));

subplot(2,1,2); hold on; grid on;
y_low  = min(y_ring_MC, [], 2);
y_high = max(y_ring_MC, [], 2);
fill([time*1e9, fliplr(time*1e9)], [y_low', fliplr(y_high')], [1 .8 .8], 'EdgeColor','none');
plot(time*1e9, mean(y_ring_MC,2), 'r', 'LineWidth', 1.5);
plot(time*1e9, y_ode45, 'k', 'LineWidth', 2);
legend('Tolerance Envelope', 'Mean MRR', 'ODE45');
xlabel('Time [ns]'); ylabel('y(t)');
title('Sensitivity to Coupling, Loss, and Detuning');

%% ---------------- Additional basic plots ----------------
% 1) Time-domain input & outputs
figure('Name','Time-Domain Signals');
subplot(3,1,1); plot(time*1e9, in_t, 'k', 'LineWidth', 2);
grid on; ylabel('Input x(t)'); xlim([0 tmax*1e9]);

subplot(3,1,2); hold on; grid on;
plot(time*1e9, y_ode45, 'k', 'LineWidth', 1.5);
plot(time*1e9, y_ring_nom, 'r--', 'LineWidth', 1.5);
legend('ODE45', 'MRR'); ylabel('Output y(t)'); xlim([0 tmax*1e9]);

subplot(3,1,3); hold on; grid on;
plot(time*1e9, abs(y_ode45).^2, 'k', 'LineWidth', 1.5);
plot(time*1e9, abs(y_ring_nom).^2, 'r--', 'LineWidth', 1.5);
legend('ODE45 |y|^2', 'MRR |y|^2'); xlabel('Time [ns]'); ylabel('|y(t)|^2'); xlim([0 tmax*1e9]);

% 2) Frequency-domain spectra
% Nominal drop-port TF for MRR chain
H_nom = ones(1,N);
for i=1:order
    beta = 2*pi*Df/(c/neff);
    H_nom = H_nom .* ((1 - r_nom(i)^2) .* alpha_nom(i) ./ (1 - r_nom(i)^2 .* alpha_nom(i) .* exp(-1j*beta*L(i))));
end
% Nominal ODE TF (ideal)
H_ode_nom = ones(1,N);
for i=1:order
    H_ode_nom = H_ode_nom .* (1/k(i)) * (1/tau_c(i)) ./ (1/tau_c(i) + 1j*2*pi*Df);
end

figure('Name','Frequency-Domain Spectra'); hold on; grid on;
plot(Df/1e9, 10*log10(abs(H_nom./max(abs(H_nom))).^2), 'r', 'LineWidth', 1.5);
plot(Df/1e9, 10*log10(abs(H_ode_nom./max(abs(H_ode_nom))).^2), 'b', 'LineWidth', 1.5);
plot(Df/1e9, 10*log10(abs(IN./max(abs(IN))).^2), 'k', 'LineWidth', 1.5);
xlabel('Frequency [GHz]'); ylabel('Magnitude [dB]');
xlim([-15 15]); ylim([-30 0]); legend('MRR TF','Ideal ODE TF','Input Spectrum');

%% ---------------- State-space ODE function ----------------
function dydt = ss_odes(~, y, k_vec, xval)
    N = numel(k_vec);
    dydt = zeros(N,1);
    dydt(1) = -k_vec(1)*y(1) + xval;
    for j = 2:N
        dydt(j) = -k_vec(j)*y(j) + y(j-1);
    end
end