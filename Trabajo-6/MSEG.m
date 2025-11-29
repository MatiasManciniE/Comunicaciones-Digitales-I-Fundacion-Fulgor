
clear 
close all
clc

h_ref = fir1(4, 0.3)'; % FILTRO DE REFERENCIA
Ncoeff=length(h_ref); % CANTIDAD DE COEFICIENTES

Lbatch=1e3; % TAMAÑO DE Xk

Niters = 10000; % Iteraciones

c = zeros(Ncoeff,1); % Respuesta inicial de C

Delta_tap = 0.0001; % Delta para estimar el gradiente
beta = 0.001; % Paso del algoritmo

% Inicializo el gradiente
gradient = zeros(Ncoeff,1);

% Variable para loging
mse_logeo = zeros(Niters,1);
gradient_logeo = zeros(Niters, Ncoeff);
c_logeo = zeros(Niters, Ncoeff);

for iter=1:Niters
    
    % Entrada de datos
    x = randn(Lbatch, 1);
    
    % Salida de referencia
    y_ref = filter(h_ref,1,x);

    % Salida de mi estimador
    y=filter(c,1,x);
    
    % Calculo el MSE
    mse0 = mean(abs(y-y_ref).^2); % Calculo la esperanza como la media en el tiempo
    
    % Loggeo el MSE para cada iteracion
    mse_logeo(iter) = mse0;
    
    % Calculo el gradiente
    for nc=1:Ncoeff % Recorro los taps del filtro
        
        % Estimacion del gradiente para el coeficiente nc
        caux = c; 
        caux(nc) = caux(nc) + Delta_tap;
        y = filter(caux,1,x);
        mse1 = mean(abs(y-y_ref).^2);
        gradient(nc)=(mse1-mse0)/Delta_tap;
        
    end
    
    % Ejecuto la formula del algoritmo
    c = c - beta*gradient;
    
    % Logeo gradiente y los taps
    gradient_logeo(iter,:)=gradient;
    c_logeo(iter,:)=c;
    
end

figure
stem(h_ref);
hold all
stem(c,'--');
title("Respuesta referencia vs respuesta estimada")
legend("Ref", "Est")

figure
plot(h_ref-c,'-o')
title("Resta referencia menos estimada")

figure
plot(c_logeo)
title("Avance de la estimación")

figure
plot(10*log10(mse_logeo))
title("Avance del MSE")

figure
plot(gradient_logeo)
title("Avance del gradiente")