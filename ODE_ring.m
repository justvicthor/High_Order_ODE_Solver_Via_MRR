%% ODE_ring.m

close all
clear all 


k = 0.5;  % ns^-1
A=1e9;  %time scaling parameter  

% Input signal x(t)
C=4;
%x = @(t) sin(t*1e12); %sinusoidal signal
%A*t*exp(-(A*t).^2).*cos(3*A*t);   %arbitrary signal
x = @(t) C*(t > 0);    % step function (Heaviside) 

% Definition of the ODE
odefun = @(t,y) A*(x(t) - k*y);

% Initial condition
y0 = 1;

% Time span
tspan = [-100e-9 100e-9];

% Solve using ode45
[t, y] = ode45(odefun, tspan, y0);


% Implementing the ODE solver with a microring resonator

k_ring=k*A;   
tau_c=1/k_ring;  %cavity life time of the RR
R=5000e-6;  %radiurn of the MRR
L_ring=2*pi*R;
c=3e8;
neff=1.5;   %effective index of the MRR waveguide
tau=L_ring/(c/neff);  %round trip time 
tau_n=tau_c/tau;
r=sqrt(tau_n/(1+tau_n));  % coupling coefficient of the directional coupler of the MRR 

%numerical version of the input step function x(t)
N=1e5;
time=linspace(min(t),max(t),N);
dt=time(2)-time(1);
in_ring = zeros(size(time));
in_ring(find(time>0))=C;
IN_ring=fftshift(fft(in_ring));

Df=linspace(-1/(2*dt),1/(2*dt),N);
beta=2*pi*Df/c*neff;
H_drop=1/k*(1-r^2)./(1-r^2*exp(-j*beta*L_ring)); %frequency domain description of the MRR
H_ODE=1/k*(1/tau_c)./(1/tau_c+j*2*pi*Df); %frequency domain descritpion of the ODE

Out_ring=IN_ring.*H_drop;
Out_ODE=IN_ring.*H_ODE;

out_ring=ifft(fftshift(Out_ring));
out_oude=ifft(fftshift(Out_ODE));

%generate plots
t_min=-1;t_max=20;

figure(1);
subplot(311);hold on;grid on;box on;plot(time*1e9, in_ring,'k','LineWidth',2)
xlabel('Time [ns]')
ylabel('Input x(t)')
xlim([t_min t_max])
set(gca,'fontsize',12)

subplot(312);hold on;grid on;box on;
plot(t*1e9, y,'k','LineWidth',2)
grid on; plot(time*1e9, out_ring,'r','LineWidth',2)
xlabel('Time [ns]')
ylabel('Output y(t)')
xlim([t_min t_max])
set(gca,'fontsize',12)

subplot(313);hold on;grid on;box on;plot(t*1e9, abs(y).^2,'k','LineWidth',2)
plot(time*1e9, (out_ring).^2,'r','LineWidth',2)
xlabel('Time [ns]')
ylabel(' Output | y(t) |^2')
xlim([t_min t_max])
set(gca,'fontsize',12)

figure(2);hold on; grid on, box on
plot(Df/1e9,10*log10(abs(H_drop./max(abs(H_drop))).^2),'r','LineWidth',2)
plot(Df/1e9,10*log10(abs(H_ODE./max(H_ODE)).^2),'b','LineWidth',2)
plot(Df/1e9,10*log10(abs(IN_ring./max(abs(IN_ring))).^2),'k','LineWidth',2)
set(gca,'fontsize',12)
ylim([-30 0])
xlim([-15 15])
xlabel('Frequency [GHz]')
ylabel('Spectrum [dB]')
%title('Solution of dy/dt + ky = x(t)')



