sim('simulacao')
t = 0:0.1:10;
N = 101;
T = 0.1;
A = [1 T;0 1];
B = [T^2/2; T];
C = [1 0];
k = [1.8; 8.1];
u = IO.signals.values(:,2);
y = IO.signals.values(:,1);
xk = [1; 0];
x1e = zeros(N,1);
x2e = zeros(N,1);
x1 = Estado.signals.values(:,1);
x2 = Estado.signals.values(:,2);
x1e(1) = xk(1);
x2e(1) = xk(2);
for i = 2:N
    xk = A*xk + B*u(i) + k*(y(i) - C*xk);
    x1e(i) = xk(1);
    x2e(i) = xk(2);
end
figure(1);
hold on;
plot(t,x1e);
plot(t,x1);

figure(2);
hold on;
plot(t,x2e);
plot(t,x2);

