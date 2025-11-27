clc; clear all; close all;

%% ===================== Parámetros generales =====================
BR = 32e9;      % Baud rate
L  = 1e5;       % Cantidad de símbolos
N  = 2;         % Oversampling
k  = 2;         % log2(M) para QPSK
M  = 4;         % QPSK

% Barrido de Eb/No (dB)
EbNo_all = 0:0.5:12;      % podés ajustarlo
nEb      = length(EbNo_all);

% Barrido de ancho de banda del canal (frecuencia de corte normalizada)
% (parámetro del fir1 respecto a fs/2)
BW_all   = [0.25 0.3 0.35 0.4 0.45];
nBW      = length(BW_all);

% Matriz para guardar BER: filas = BW, columnas = Eb/No
BER_A = zeros(nBW, nEb);

%% ===================== Barrido de simulaciones =====================
for iBW = 1:nBW
    BWnorm = BW_all(iBW);
    fprintf('Simulando BWnorm = %.2f (caso A: BW antes del ruido)\n', BWnorm);
    
    for iEb = 1:nEb
        EbNo = EbNo_all(iEb);
        
        BER_A(iBW, iEb) = sim_FSE_caseA(BR, L, N, k, M, EbNo, BWnorm);
        
        fprintf('   EbNo = %2d dB  ->  BER = %e\n', EbNo, BER_A(iBW,iEb));
    end
end

%% ===================== Plots BER vs EbNo (caso A) =====================
figure;
for iBW = 1:nBW
    semilogy(EbNo_all, BER_A(iBW,:), 'o-', 'LineWidth', 1.3); hold on;
end
grid on;
xlabel('E_b/N_0 [dB]');
ylabel('BER');
legend(cellstr(num2str(BW_all','BW = %.2f')));
title('Caso A: Limitación de BW ANTES del ruido (sistema FSE)');

%% ===================== Penalidad de SNR vs BW =====================
% Definimos una BER objetivo
BER_target = 1e-3;

% Eb/No teórico de referencia para QPSK en AWGN
EbNo_lin_ref = (qfuncinv(BER_target)).^2 / 2;   % BPSK/QPSK
EbNo_dB_ref  = 10*log10(EbNo_lin_ref);

penalty_A = nan(nBW,1);

for iBW = 1:nBW
    ber_curve = BER_A(iBW,:);
    
    if any(ber_curve <= BER_target)
        % Interpolamos en escala log(BER) para encontrar EbNo requerido
        Eb_req_A = interp1(log10(ber_curve), EbNo_all, ...
                           log10(BER_target), 'linear');
        penalty_A(iBW) = Eb_req_A - EbNo_dB_ref;
    end
end

figure;
plot(BW_all, penalty_A, 'o-', 'LineWidth', 1.5); grid on;
xlabel('Ancho de banda del canal (W_c normalizado)');
ylabel('Penalidad de SNR [dB]');
title(sprintf('Caso A: Penalidad de SNR para BER = %.1e', BER_target));

function ber = sim_FSE_caseA(BR, L, N, k, M, EbNo, BWnorm)

    % --------- PARÁMETROS DERIVADOS ---------
    fs   = N*BR;
    BETA = 0.2;
    SPAN = 25;
    SPS  = fs/BR;

    % --------- GENERACIÓN DE SÍMBOLOS QPSK ---------
    xi = 2*randi([0,1],L,1)-1;
    xq = 2*randi([0,1],L,1)-1;
    x  = xi + 1j.*xq;             % símbolos QPSK a tasa de símbolo

    % Upsampling
    xup = upsample(x,N);

    % Filtro de Tx (RRC)
    g = rcosdesign(BETA, SPAN, SPS, "sqrt");
    g = g/sum(g);
    yup = filter(g,1,xup);

    % --------- CANAL con limitación de BW ANTES del ruido ---------
    b   = fir1(10, BWnorm);       % ancho de banda del canal (normalizado)
    ych = filter(b,1,yup);        % filtrado por el canal

    % --------- RUIDO AWGN ---------
    SNR_dB = EbNo - 10*log10(N) + 10*log10(k);
    SNR    = 10^(SNR_dB/10);
    Pt     = var(yup);
    No     = Pt/SNR;
    sigma  = sqrt(No/2);
    ruido  = sigma.*(randn(length(yup),1) + 1j.*randn(length(yup),1));

    rx = ych + ruido;             % *** BW antes del ruido (caso A) ***

    % --------- ECUALIZADOR FSE (LMS) ---------
    NTAPS = 101;
    step  = 2e-3;
    leak  = 1e-3;

    Xbuffer      = zeros(NTAPS,1);
    W            = zeros(1,NTAPS);
    W((NTAPS+1)/2) = 1.0;

    LY          = length(rx);
    out_eq_down = zeros(floor(LY/N),1);

    idx_sym = 0;
    for m = 1:LY-NTAPS-1
        Xbuffer(2:end) = Xbuffer(1:end-1);
        Xbuffer(1)     = rx(m);
        
        yeq = W*Xbuffer;
        
        if mod(m,N)==0
            idx_sym = idx_sym + 1;
            out_eq_down(idx_sym) = yeq;
            
            errorLMS = yeq - slicer_qam(yeq,M);
            grad     = errorLMS .* Xbuffer';
            W        = W*(1-step*leak) - step*grad;
        end
    end

    % --------- CÁLCULO DE BER con my_ber_checker ---------
    ak_hat_v = out_eq_down(out_eq_down ~= 0);
    ak_hat_v = ak_hat_v(:);

    Ltx  = length(x);
    Lrx  = length(ak_hat_v);
    Lmin = min(Ltx, Lrx);

    ref_v  = x(1:Lmin);         % símbolos TX
    data_v = ak_hat_v(1:Lmin);  % símbolos RX ecualizados

    [ber, ~] = my_ber_checker(data_v, ref_v, M, 'auto');
end
