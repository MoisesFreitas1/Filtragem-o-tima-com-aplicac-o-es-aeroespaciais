clc
clear all
close all

x_ = [0;0];
P_ = eye(2);
T = 100;
Ts = 0.1;
Nr = 1;
Q = 0.01*eye(2);
R = 0.01;
Hk1 = [1 0];


xe1kk = zeros(2,T/Ts);
xtil = zeros(2,T/Ts);
ye1kk = zeros(1,T/Ts);
P1kk = zeros(2,2,T/Ts);
Pkky = zeros(1,T/Ts);
Pkkxy = zeros(2,T/Ts);
K1k = zeros(2,T/Ts);
sigma1 = zeros(1,T/Ts); % Desvio padrão de X1
sigma2 = zeros(1,T/Ts); % Desvio padrão de X2

xe1kk(:,1) = x_;
P1kk(:,:,1) = P_;

for i = 1:Nr
    x0 = x_ + sqrt(P_)*randn(2,1);
    sim('simulacaoEKF');
    y1k = Dadosy.signals.values';
    uk = Dadosu.signals.values';
    x1k = Dadosx.signals.values';
    for j = 1:T/Ts
        %Predição
        [xe1kk(:,j),P1kk(:,:,j)] = integrador(xe1kk(:,j),P1kk(:,:,j),uk(j),Q,Ts);
        Pkky(:,j) = Hk1*P1kk(:,:,j)*Hk1' + R;
        Pkkxy(:,j) = P1kk(:,:,j)*Hk1';
        hk1 = xe1kk(1,j);
        
        %Atualização
        ye1kk(:,j) = hk1;
        K1k(:,j) = Pkkxy(:,j)*inv(Pkky(:,j));
        xe1kk(:,j+1) = xe1kk(:,j) + K1k(:,j)*(y1k(:,j) - ye1kk(:,j));
        P1kk(:,:,j+1) = P1kk(:,:,j) - Pkkxy(:,j)*inv(Pkky(:,j))*Pkkxy(:,j)';
        
        %Erro
        xtil(:,j) = x1k(:,j) - xe1kk(:,j);
        
        %Desvio-padrão
        sigma1(1,j) = sqrt(P1kk(1,1,j));
        sigma2(1,j) = sqrt(P1kk(2,2,j));
    end
    
    figure(1);
    hold on;
    plot(xe1kk(1,:),'r');
    plot(x1k(1,:),'b');
    title('X1');
    
    figure(2);
    hold on;
    plot(xe1kk(2,:),'r');
    plot(x1k(2,:),'b');
    title('X2');
    
    figure(3);
    hold on;
    plot(xtil(1,:),'b');
    plot(sigma1(1,:),'r');
    plot(-sigma1(1,:),'r');
    
    figure(4);
    hold on;
    plot(xtil(2,:),'b');
    plot(sigma2(1,:),'r');
    plot(-sigma2(1,:),'r');
end

