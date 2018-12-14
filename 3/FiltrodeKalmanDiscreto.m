%Modelo
A = [1 0.1;0 1];
B = [0.005;0.1];
C = [1 0];
Q = [0.01 0;0 0.04];
R = 0.01;
P_ = [1e-4 0;0 1e-8];
x_ = [1;0];
yk_0 = 5;
Ts = 0.1;
T = 30;
Nr = 1;
e1 = [1;0];
e2 = [0;1];

x1k = zeros(2,T/Ts); %Estado verdadeiro do sistema
y1k = zeros(1,T/Ts); %Saída verdadeira do sistema

xe1kk = zeros(2,T/Ts); %Estado estimado pelo Filtro de Kalman
ye1kk = zeros(1,T/Ts); %Saída estimada pelo Filtro de Kalman
P1kk = zeros(2,2,T/Ts); %Variância
Pkky = zeros(1,T/Ts); %Variância da saída
Pkkxy = zeros(2,T/Ts); %Variância cruzada
K1k = zeros(2,T/Ts); %Ganho de Kalman
erro = zeros(2,T/Ts); %Erros de X1 e X2
sigma1 = zeros(1,T/Ts); % Desvio padrão de X1
sigma2 = zeros(1,T/Ts); % Desvio padrão de X2

for j = 1:Nr %Monte Carlo
    %Simulação
    x1k(:,1) = x_ + sqrt(P_)*randn(2,1);
    vk = sqrt(R)*randn;
    y1k(1) = C*x1k(:,1) + vk;
    
    %Filtro de Kalman - Inicialização
    xe1kk(:,1) = x_;
    P1kk(:,:,1) = P_;
    ye1kk(:,1) = C*xe1kk(:,1);
    P1kk(:,:,1) = A*P1kk(:,:,1)*A'+Q;
    Pkky(:,1) = C*P1kk(:,:,1)*C'+R;
    Pkkxy(:,1) = P1kk(:,:,1)*C';
    K1k(:,1) = Pkkxy(:,1)*inv(Pkky(:,1));
    xe1kk(:,1) = xe1kk(:,1) + K1k(:,1)*(y1k(1)-ye1kk(:,1));
    P1kk(:,:,1) = P1kk(:,:,1) - Pkkxy(:,1)*inv(Pkky(:,1))*Pkkxy(:,1)';
    %Cálculo do erro
    erro(1,1) = e1'*(x1k(:,1)-xe1kk(:,1)); %erro de X1
    erro(2,1) = e2'*(x1k(:,1)-xe1kk(:,1)); %erro de X2
    %Desvio padrão
    sigma1(1,1) = sqrt(P1kk(1,1,1));
    sigma2(1,1) = sqrt(P1kk(2,2,1));
    for i = 2:T/Ts
        %Calculo do x1k
        wk = sqrt(Q)*randn(2,1);
        uk = 10*(yk_0 - e1'*x1k(:,i-1)) - 2*e2'*x1k(:,i-1);
        x1k(:,i) = A*x1k(:,i-1) + B*uk + wk;
        %Calculo do y1k
        vk = sqrt(R)*randn;
        y1k(i) = C*x1k(:,i) + vk;
        
        %Filtro de Kalman - Predição
        xe1kk(:,i) = A*xe1kk(:,i-1) + B*uk;
        ye1kk(:,i) = C*xe1kk(:,i);
        P1kk(:,:,i) = A*P1kk(:,:,i-1)*A'+Q;
        Pkky(:,i) = C*P1kk(:,:,i)*C'+R;
        Pkkxy(:,i) = P1kk(:,:,i)*C';
        
        %Filtro de Kalman - Atualização
        K1k (:,i) = Pkkxy(:,i)*inv(Pkky(:,i));
        xe1kk(:,i) = xe1kk(:,i) + K1k(:,i)*(y1k(i)-ye1kk(:,i));
        P1kk(:,:,i) = P1kk(:,:,i) - Pkkxy(:,i)*inv(Pkky(:,i))*Pkkxy(:,i)';
        
        %Cálculo do erro
        erro(1,i) = e1'*(x1k(:,i)-xe1kk(:,i)); %erro de X1
        erro(2,i) = e2'*(x1k(:,i)-xe1kk(:,i)); %erro de X2
        %Desvio padrão
        sigma1(1,i) = sqrt(P1kk(1,1,i));
        sigma2(1,i) = sqrt(P1kk(2,2,i));
    end
        
    yk_ = zeros(1,T/Ts);
    for i = 1:T/Ts
        yk_ (1,i) = 5;
    end
    
    %Posição vertical
    figure(1);
    hold on;
    plot(y1k);
    plot(ye1kk);
    plot(yk_);
    title('Posição vertical');

    %X1
    figure(2);
    hold on;
    plot(xe1kk(1,:));
    plot(x1k(1,:));
    title('x1');
    
    %X2
    figure(3);
    hold on;
    plot(xe1kk(2,:));
    plot(x1k(2,:));
    title('x2');
    
    %Erro de X1
    figure(4);
    hold on;
    plot(erro(1,:));
    plot(sigma1(1,:));
    plot(-sigma1(1,:));
    title('Erro de x1');

    %Erro de X2
    figure(5);
    hold on;
    plot(erro(2,:));
    plot(sigma2(1,:));
    plot(-sigma2(1,:));
    title('Erro de x2'); 
end