clc;
clear;

teta_true = [1;2]; %teta verdadeiro
R = 0.01*[1 0;0 1]; %matriz de covariância
R_inv = inv(R); %inversa da matriz de covariância
Nr = 10000; %número de realizações
N = 50; %número de somas em teta
m = zeros(2,1); %media que inicializa zerada
e = zeros(2,1); %erro que inicializa zerado
RMS = zeros(2,1); %valor RMS que inicializa zerado
teta_s = zeros(2,Nr); %variável para armazenar todos os valores de teta estimado
for i = 1:Nr
    S1 = zeros(2,2); %soma para calcular teta. Essa soma não depende da saída do sistema. Devem ser zeradas a cada vez que se calcula um teta, por esse motivo está dentro deste for
    S2 = zeros(2,1); %soma para calcular teta. Devem ser zeradas a cada vez que se calcula um teta, por esse motivo está dentro deste for
    for j = 1:N
        h = [1 j;j sin(j/200*pi)];
        v = sqrt(R)*randn(2,1);
        y = h*teta_true+v; %cálculo da saída do sistema. Imagine que você tem um sensor capaz de medir esse parâmetro. Ele está sendo simulado aqui.
        
        S1 = S1 + h'*R_inv*h; %soma independente das medições
        S2 = S2 + h'*R_inv*y; %soma dependente das medições
    end
    teta_e = inv(S1)*S2; %Calcula um teta estimado
    teta_s(1,i) = teta_e(1); %armazena a primeira componente de teta
    teta_s(2,i) = teta_e(2); %armazana a segunda componente de teta

    m = m + teta_e;
    e = teta_true - teta_e; %erro
    e = e.^2; %erro quadrático
    RMS = RMS + e;
end
RMS_te = sqrt(inv(S1)) %RMS teórico
m = m/Nr %média
RMS = sqrt(RMS/Nr) %RMS experimental

figure(1);
histogram(teta_s(1,:)) %plota o histograma da primeira componente de teta
hold on
plot(m(1),0,'r*') %plota a média da primeira componente de teta
plot(m(1)-RMS(1),0,'g*') %plota a média menos o RMS da primeira componente de teta
plot(m(1)+RMS(1),0,'g*') %plota a media mais o RMS da primeira componente de teta

figure(2);
histogram(teta_s(2,:)) %plota o histograma da segunda componente de teta
hold on
plot(m(2),0,'r*') %plota a média da segunda componente de teta
plot(m(2)-RMS(2),0,'g*') %plota a média menos o RMS da segunda componente de teta
plot(m(2)+RMS(2),0,'g*') %plota a media mais o RMS da segunda componente de teta