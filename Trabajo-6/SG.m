clear all
close all

rand('seed',1)
h_ref = fir1(4,0.3)';
Ncoeff=length(h_ref);
Lsim=10e3;

x = randn(Lsim, 1);
y_ref = filter(h_ref,1,x) + 0.0025*randn(Lsim,1);

c=zeros(Ncoeff,1);
step = 1e-1;
buffer=zeros(size(c));

y=zeros(size(x));
error=zeros(size(x));
mse_log=zeros(Lsim,1);
c_log=zeros(Lsim, Ncoeff);

for n=1:Lsim
    
    buffer(2:end)=buffer(1:end-1);
    buffer(1)=x(n);
    y(n) = sum(c.*buffer);
    
    error(n) = y(n)-y_ref(n);
    c = c - step*error(n)*buffer;
        
    c_log(n,:)=c;
    
end

figure;
plot(10*log10(filter(ones(2560,1)/2560,1,abs(error).^2)))

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
plot(c_log)
title("Avance de la estimación")