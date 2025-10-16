%-----------------------------------------------------------------------------%
%                               Ejercicio 8 (Corregido y Mejorado)
%-----------------------------------------------------------------------------%
clear; close all; clc;

%% Parámetros básicos
L = 10000;          % cantidad de símbolos
BR = 16e9;          % baud rate
N  = 4;             % sobremuestreo
rolloff = 0.5;      % roll-off
span = 6;           % span de símbolos del filtro
fs = N*BR;          % frecuencia de muestreo
Ts = 1/fs;

%% 1. Generar señales I y Q independientes (PAM2)
I = 2*randi([0,1],L,1)-1;
Q = 2*randi([0,1],L,1)-1;
m = I + 1j*Q;       % señal compleja QPSK

% Sobremuestreo
m_up = upsample(m,N);

%% 2. Filtrado (root raised cosine)
h = rcosdesign(rolloff, span, N, 'sqrt');
m_filt = conv(m_up,h,'same');

%% 3. Modulador analítico
f0 = (1+rolloff)*BR/2 * 1.2;
t = (0:length(m_filt)-1)'*Ts;
s = m_filt .* exp(1j*2*pi*f0*t);

%% 4. Demodulación ideal
g = s .* exp(-1j*2*pi*f0*t);

%% 5. Compensación de retardo y decimación
delay = span*N/2;                   % retardo de grupo
g_sync = g(delay+1:end);
m_sync = m_filt(delay+1:end);

% Decimación (1 muestra por símbolo)
g_dec = g_sync(1:N:end);
m_dec = m_sync(1:N:end);

%% a) Constelación transmitida
figure('Color','k','Name','Constelación Transmitida');
scatter(real(m_dec), imag(m_dec), 10, 'y', 'filled');
title('Constelación transmitida m(t)', 'Color','w');
xlabel('In-Phase','Color','w'); ylabel('Quadrature','Color','w');
grid on;

%% b) Constelación recibida sin error
figure('Color','k','Name','Constelación Recibida');
scatter(real(g_dec), imag(g_dec), 10, 'y', 'filled');
title('Constelación recibida g(t) - sin error', 'Color','w');
xlabel('In-Phase','Color','w'); ylabel('Quadrature','Color','w');
grid on;

%% c) Error de fase constante
theta_list = [pi/8, pi/4, pi/2];
for k = 1:length(theta_list)
    theta = theta_list(k);
    p = exp(-1j*theta);
    g_phase = g .* p;
    g_dec_phase = g_phase(delay+1:end);
    g_dec_phase = g_dec_phase(1:N:end);

    figure('Color','k','Name',['Error de fase ',num2str(theta/pi),'π']);
    scatter(real(g_dec_phase), imag(g_dec_phase), 10, 'y', 'filled');
    title(['Constelación recibida con error de fase θ = ',num2str(theta/pi),'·π rad'], 'Color','w');
    xlabel('In-Phase','Color','w'); ylabel('Quadrature','Color','w');
    grid on;
end

%% d) Comparar θ = 0 vs θ = π/2
theta = [0, pi/2];
colors = ['b','r'];
figure('Color','k','Name','Comparación de fases'); hold on; grid on;
for k=1:2
    p = exp(-1j*theta(k));
    g_phase = g .* p;
    g_dec_phase = g_phase(delay+1:end);
    g_dec_phase = g_dec_phase(1:N:end);
    scatter(real(g_dec_phase), imag(g_dec_phase), 10, colors(k), 'filled');
end
legend('θ = 0','θ = π/2','TextColor','w','Location','best');
title('Comparación de constelaciones con diferentes fases', 'Color','w');
xlabel('I','Color','w'); ylabel('Q','Color','w');

%% e) Offset de frecuencia Δf
df = 0.02*BR;  % 2% del baud rate (~320 MHz para BR=16GHz)
p_freq = exp(-1j*2*pi*df*t);
g_off = g .* p_freq;
g_off_dec = g_off(delay+1:end);
g_off_dec = g_off_dec(1:N:end);

figure('Color','k','Name','Offset de frecuencia');
scatter(real(g_off_dec), imag(g_off_dec), 10, 'y', 'filled');
title(['Constelación con offset de frecuencia Δf = ', num2str(df/1e6), ' MHz'], 'Color','w');
xlabel('In-Phase','Color','w'); ylabel('Quadrature','Color','w');
grid on;

%-----------------------------------------------------------------------------%
% Comentarios:
% - (a) La constelación transmitida
