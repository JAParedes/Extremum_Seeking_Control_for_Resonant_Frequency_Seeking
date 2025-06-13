% Program for plotting the magnitude of the frequency response of the
% mass-spring-damper system using the analytical expression derived in the
% attached document 

close all
clear all
clc

%% Setup

m = 1;    % Mass
k = 1000;    % Spring stiffness coefficient
c = 2;      % Damper dampening coefficient
b0 = 1/m;   % Input gain

omega_n = sqrt(k/m);                    % Natural frequency of mass-spring-damper system
zeta = c/(2*sqrt(m*k));                 % Damping ratio of mass-spring-damper system
omega_r = omega_n*sqrt(1-2*zeta^2);     % Resonant frequency of mass-spring-damper system

ww = 0:0.01:60; % Frequency range for magnitude evaluation

% Magnitude of mass-spring-damper system transfer function using the
% analytical solution
magG = abs(b0)./sqrt(ww.^4 - 2*(omega_n^2)*(1 - 2*zeta^2).*(ww.^2) + omega_n^4);

%% Plotting

fontLabels = 24;

figure(1)

set(gcf, 'color', [1 1 1])
set(gcf, 'position', [680, 206, 850, 450])

plot(ww,magG, 'LineWidth', 2)
if isreal(omega_r) 
    xline(omega_r, 'r--', 'LineWidth', 2)
end

set(gca,'TickLabelInterpreter','latex')
ax = gca;
ax.FontSize = fontLabels;

ylabel('Magnitude', 'interpreter', 'latex', 'fontsize', fontLabels)
xlabel('$\omega$ (rad/s)', 'interpreter', 'latex', 'fontsize', fontLabels)

grid on
box on

%% Comparison with magnitude obtained from 'bode' command.

G_tf = tf(b0,[1 2*zeta*omega_n omega_n^2]);
[magG2, ~, ~] = bode(G_tf, ww);

figure(2)

set(gcf, 'color', [1 1 1])
set(gcf, 'position', [680, 206, 850, 450])

plot(ww,magG, 'LineWidth', 2)
hold on
plot(ww,squeeze(magG2),'--', 'LineWidth', 2)
hold off
if isreal(omega_r) 
    xline(omega_r, 'r--', 'LineWidth', 2)
end

legend({'Analytical Solution','Solution using ``bode" command'},'interpreter','latex','fontsize', fontLabels);

set(gca,'TickLabelInterpreter','latex')
ax = gca;
ax.FontSize = fontLabels;

ylabel('Magnitude', 'interpreter', 'latex', 'fontsize', fontLabels)
xlabel('$\omega$ (rad/s)', 'interpreter', 'latex', 'fontsize', fontLabels)

grid on
box on