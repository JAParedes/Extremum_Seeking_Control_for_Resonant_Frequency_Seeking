% Program tests the ability of the proposed nonlinear oscillators to switch
% frequencies

close all
clear all
clc

%% Setup

% First frequency
f1 = 20; %Hz
omega1 = 2*pi*f1;

% Second frequency
f2 = 50; %Hz
omega2 = 2*pi*f2;

amp = 1;        % Output amplitude.
lambda = 1;     % Amplitude convergence rate factor.

Ts = 0.001;     % Sampling time.
Tfin = 0.5;     % Final simulation time.
t_switch = 0.25;% Time for frequency switch.

x0 = [0.0; 1.0];% Initial state

%% Simulations

% Continous-time nonlinear oscillator simulation
[tt,xx_CT] = ode45(@(t,x)nonl_osc_CT(t, x, step_fun(t, t_switch, omega1, omega2), amp, lambda), 0:Ts:Tfin, x0);

tt = tt.';
xx_CT = xx_CT.';
xx_DT = zeros(size(xx_CT));
xx_DT(:,1) = x0;

% Discrete-time nonlinear oscillator simulation
for ii = 2:length(tt)
    omega = step_fun(Ts*(ii-1), t_switch, omega1, omega2);
    xx_DT(:,ii) = nonl_osc_DT(xx_DT(:,ii-1), omega, amp, Ts, lambda);
end

%% Plotting

fontLabels = 18;

figure(1)

set(gcf, 'color', [1 1 1])
set(gcf, 'position', [680, 206, 850, 450])

plot(tt,xx_CT(1,:), 'LineWidth', 2)
hold on
stairs(tt,xx_DT(1,:), 'LineWidth', 2)
hold off

legend({'Continuous-time Nonlinear Occillator','Discrete-time Nonlinear Occillator'},'interpreter','latex','fontsize', fontLabels);

set(gca,'TickLabelInterpreter','latex')
ax = gca;
ax.FontSize = fontLabels;

ylabel('$y$', 'interpreter', 'latex', 'fontsize', fontLabels)
xlabel('$t$ (s)', 'interpreter', 'latex', 'fontsize', fontLabels)

grid on
box on

% Continous-time nonlinear oscillator dynamics
function xdot = nonl_osc_CT(t, x, omega, amp, lambda)
    xdot = ([0 omega; -omega 0] - lambda*((x.'*x)- amp^2)*eye(2))*x;
end

% Discrete-time nonlinear oscillator update equation
function xkp1 = nonl_osc_DT(xk, omega, amp, Ts, lambda)
    xkp1 = ([cos(omega*Ts) sin(omega*Ts); -sin(omega*Ts) cos(omega*Ts)] - lambda*((xk.'*xk)- amp^2)*eye(2))*xk;
end

% Step function for frequency switching.
function out = step_fun(t, t_switch, val1, val2)
    out = double(t<=t_switch).*val1 + double(t>t_switch).*val2;
end