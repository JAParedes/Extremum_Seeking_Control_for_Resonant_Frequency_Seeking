% Program for running extremum seeking for resonant frequency seeking in
% discrete-time. This simulation example is applied to the
% mass-spring-damper system. In this example, the demodulation dither
% signal is set to zero after reaching a neighborhood near the resonant
% frequency.

close all
clear all
clc

%% Setup

Tfin = 400;         % Final simulation time
Ts_sim = 0.001;     % Simulation sampling time
%Ts_sim = 0.0001;
Ts = 0.001;         % Sampling time for nonlinear oscillator and moving RMS filter
Ts_esc = 0.1;       % Sampling time for DT/ESC

% Mass-spring-damper system parameters
m = 1;      % Mass
k = 1000;   % Spring stiffness coefficient
c = 2;      % Damper damping coefficient

x0 = [0;0]; % Initial mass-spring-damper system state

omega_n = sqrt(k/m);                % Natural frequency of mass-spring-damper system
zeta = c/(2*sqrt(m*k));             % Damping ratio of mass-spring-damper system
omega_r = omega_n*sqrt(1-2*zeta^2); % Resonant frequency of mass-spring-damper system

% Discrete-time, nonlinear oscillator parameters
A_cmd = 1;              % Oscillator commanded amplitude
lambda_osc = 1;         % Amplitude convergence rate factor
x_osc_0 = [0; A_cmd];   % Initial oscillator state
omega_cmd_0 = 20;       % Initial oscillator commanded frequency

% Moving root mean square (RMS) parameters
l_rms = 2000;   % Number of samples for the cumulative averaging filter
x_rms_0 = 0.00; % Initial RMS filter state

% Discrete-time extremum seeking control hyperparameters
K_g = 500;      % Measurement gain
omega_l = 4;   % Low-pass filter cutoff frequency in rad/s
omega_h = 0.1;    % High-pass filter cutoff frequency in rad/s
K_esc = 2;      % Modulation gain
omega_esc = 1;  % Dither signal frequency in rad/s
A_esc_m = 1.2;  % Modulation dither signal amplitude
A_esc_d = 1.2;  % Demodulation dither signal initial amplitude

t_en = 2*(l_rms + 1)*Ts;        % Time at which ESC is enabled

t_stop = 300;                   % Time after which the modulation dither becomes zero
k_stop = round(t_stop/Ts_esc);  % ESC step after which the modulation dither becomes zero

%% Running the simulation in Simulink

out = sim('msd_resonant_frequency_seeking_DT.slx');

tt = out.tout;
force = out.force;
y = out.y;
y_rms = out.y_rms;
omega_cmd = out.omega_cmd;

%% Plotting

fontLabels = 16;

figure(1)

set(gcf, 'color', [1 1 1]) 
set(gcf, 'position', [300,300,850,420])

pp1 = subplot(3,1,1);
plot(tt,y)
ax = gca;
ax.FontSize = fontLabels; 
set(gca,'TickLabelInterpreter','latex')

xlim([0 Tfin])

grid on
box on
ylabel('$y$ (m)','interpreter','latex','fontsize',fontLabels)
xticklabels('')

pp2 = subplot(3,1,2);
plot(tt,y_rms)
ax = gca;
ax.FontSize = fontLabels; 
set(gca,'TickLabelInterpreter','latex')

xlim([0 Tfin])
ylim([0 0.012])

grid on
box on
ylabel('$y_{\rm rms}$ (m)','interpreter','latex','fontsize',fontLabels)
xticklabels('')

pp3 = subplot(3,1,3);
yline(omega_r,'r--','LineWidth',2)
hold on
plot(tt,omega_cmd)
hold off
ax = gca;
ax.FontSize = fontLabels; 
set(gca,'TickLabelInterpreter','latex')

xlim([0 Tfin])

grid on
box on
ylabel('$\omega_{\rm cmd}$ (rad/s)','interpreter','latex','fontsize',fontLabels)
xlabel('$t$ (s)','interpreter','latex','fontsize',fontLabels)

%%
set(pp1,'position',[0.11, 0.74, 0.85, 0.23])
set(pp2,'position',[0.11, 0.44, 0.85, 0.23])
set(pp3,'position',[0.11, 0.14, 0.85, 0.23])

%%
ww = 10:0.01:40; % Frequency range for magnitude evaluation

% Magnitude of mass-spring-damper system transfer function using the
% analytical solution
magG = abs(1/m)./sqrt(ww.^4 - 2*(omega_n^2)*(1 - 2*zeta^2).*(ww.^2) + omega_n^4);
% Magnitude of the closed-loop simulation at each time step
mag_CT_ESC = abs(1/m)./sqrt(omega_cmd.^4 - 2*(omega_n^2)*(1 - 2*zeta^2).*(omega_cmd.^2) + omega_n^4);

c_color = linspace(0,1,length(mag_CT_ESC));

figure(2)

set(gcf, 'color', [1 1 1])
set(gcf, 'position', [680,460,850,360])

plot(ww,magG, 'LineWidth', 2)
hold on
scatter(omega_cmd,mag_CT_ESC, 30, c_color,'filled')
hold off
xline(omega_r, '--', 'Color', [0 0.5 0],'LineWidth', 2)

set(gca,'TickLabelInterpreter','latex')
ax = gca;
ax.FontSize = fontLabels;

ylabel('Magnitude', 'interpreter', 'latex', 'fontsize', fontLabels)
xlabel('$\omega$ (rad/s)', 'interpreter', 'latex', 'fontsize', fontLabels)

grid on
box on

colormap('jet')
cb = colorbar('Ticks',[0, 1],'TickLabels',{'0','$T_{\rm fin}$'},'TickLabelInterpreter','latex','FontSize',fontLabels);
cb.Label.String = '$t$ (s)';
cb.Label.Interpreter = 'latex';