%-----------------------------------------------------------------------------%
%                               Ejercicio 7 - QPSK / QAM-4
%-----------------------------------------------------------------------------%
clear; close all; clc;

%% Parámetros básicos
L = 10000;       % número de símbolos
BR = 16e9;       % baud rate
N  = 4;          % sobremuestreo
rolloff = 0.5;   % roll-off del raised cosine
fs = N*BR;       % frecuencia de muestreo
Ts = 1/fs;

%% 1. Generación de dos señales PAM2 independientes
I = 2*randi([0,1],L,1)-1;   % Rama en fase
Q = 2*randi([0,1],L,1)-1;   % Rama en cuadratura

% Señal compleja m(t)
m = I + 1j*Q;

% Sobremuestreo
m_up = upsample(m, N);

%% 2. Filtro raised-cosine (pulso de forma)
h = rcosdesign(rolloff, 6, N, 'sqrt');
m_filt = conv(m_up, h, 'same');

%% 3. Graficar PSD de Re{m(t)}, Im{m(t)} y m(t)
NFFT = 2048;

[PSD_I,f_I] = pwelch(real(m_filt), hanning(NFFT/2), 0, NFFT, fs);
[PSD_Q,f_Q] = pwelch(imag(m_filt), hanning(NFFT/2), 0, NFFT, fs);
[PSD_M,f_M] = pwelch(m_filt,        hanning(NFFT/2), 0, NFFT, fs);

figure('Color','w','Name','PSD de las componentes')
plot(f_I,10*log10(PSD_I/max(PSD_I)),'b','LineWidth',1.3); hold on;
plot(f_Q,10*log10(PSD_Q/max(PSD_Q)),'r','LineWidth',1.3);
plot(f_M,10*log10(PSD_M/max(PSD_M)),'k','LineWidth',1.3);
grid on
xlabel('Frecuencia [Hz]'); ylabel('PSD [dB/Hz]');
title('PSD de Re\{m(t)\}, Im\{m(t)\} y m(t)');
legend('Re\{m(t)\}','Im\{m(t)\}','|m(t)|');
xlim([-40e9 40e9]);

%% 4. Señal analítica s(t) = m(t)*e^{j2πf0t}
f0 = (1+rolloff)*BR/2 * 1.2;  % desplazamiento
t = (0:length(m_filt)-1)'*Ts;
s = m_filt .* exp(1j*2*pi*f0*t);

% PSD de s(t)
[PSD_S,f_S] = pwelch(s, hanning(NFFT/2), 0, NFFT, fs);

figure('Color','w','Name','PSD s(t)')
plot(f_S,10*log10(PSD_S/max(PSD_S)),'m','LineWidth',1.3);
xlabel('Frecuencia [Hz]'); ylabel('PSD [dB/Hz]');
title('PSD de la señal analítica s(t)');
grid on; xlim([-50e9 50e9]);

%% 5. Recuperación g(t) = s(t)*e^{-j2πf0t}
g = s .* exp(-1j*2*pi*f0*t);

% Comparar parte real e imaginaria con m(t)
figure('Color','w','Name','Comparación g(t) vs m(t)')
subplot(2,1,1);
plot(real(m_filt(1:500)),'b','LineWidth',1.3); hold on;
plot(real(g(1:500)),'r--','LineWidth',1.3);
legend('Re\{m(t)\}','Re\{g(t)\}'); grid on;
ylabel('Amplitud'); title('Parte Real');

subplot(2,1,2);
plot(imag(m_filt(1:500)),'b','LineWidth',1.3); hold on;
plot(imag(g(1:500)),'r--','LineWidth',1.3);
legend('Im\{m(t)\}','Im\{g(t)\}'); grid on;
xlabel('Muestras'); ylabel('Amplitud'); title('Parte Imaginaria');

%% Comentarios finales:
% - Re{m(t)} e Im{m(t)} tienen igual PSD y mismo ancho de banda: (1+rolloff)*BR/2 = 12 GHz.
% - m(t) ocupa el mismo ancho de banda total.
% - s(t) es analítica: sólo tiene banda positiva.
% - Aunque transmite el doble de información (I y Q), el ancho de banda se mantiene
%   porque ambas ramas se transmiten en cuadratura (ortogonales).
% - g(t) = m(t): la recuperación es exacta para partes real e imaginaria.
%-----------------------------------------------------------------------------%
