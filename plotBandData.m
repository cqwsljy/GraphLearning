% clc
% N = 500;
% theta1 = 2*pi*rand(N,1);
% band = 1.5;
% r1 = 8 + band*randn(N,1);
% r2 = 20 + band*randn(N,1);
% r3 = 30 + band*randn(N,1);
% y = [r1.*sin(theta1) ;r2.*sin(theta1);r3.*sin(theta1)];
% x = [r1.*cos(theta1) ;r2.*cos(theta1);r3.*cos(theta1)];
% 
% h = 5500;
% Knears = 10;
% X = [x,y];
% [L,d,lambda_max]=GenerateGraph_fun(X',h,Knears,'ZM'); 
% L = full(L);
%  画圆环图

[u,v] = eig(L);
[r,c] = sort(diag(v));
x2 = u(:,1);
y2 = u(:,3);

subplot(1,3,1)
plot(x(1:N),y(1:N),'*'); hold on
plot(x(N+1:2*N),y(N+1:2*N),'*');
plot(x(2*N+1:end),y(2*N+1:end),'*');

subplot(1,3,3)
plot(x2(1:N),y2(1:N),'*'); hold on
plot(x2(N+1:2*N),y2(N+1:2*N),'*');
plot(x2(2*N+1:end),y2(2*N+1:end),'*');

subplot(1,3,2)
plot(1:N,y2(1:N),'*'); hold on
plot(N+1:2*N,y2(N+1:2*N),'*');
plot(2*N+1:3*N,y2(2*N+1:end),'*');
axis([1,3*N -0.05 0.05])

figure
subplot(1,3,1)
y2 = u(:,1);
plot(1:N,y2(1:N),'*'); hold on
plot(N+1:2*N,y2(N+1:2*N),'*');
plot(2*N+1:3*N,y2(2*N+1:end),'*');
axis([1,3*N -0.05 0.05])
subplot(1,3,2)
y2 = u(:,2);
plot(1:N,y2(1:N),'*'); hold on
plot(N+1:2*N,y2(N+1:2*N),'*');
plot(2*N+1:3*N,y2(2*N+1:end),'*');
axis([1,3*N -0.05 0.05])
subplot(1,3,3)
y2 = u(:,3);
plot(1:N,y2(1:N),'*'); hold on
plot(N+1:2*N,y2(N+1:2*N),'*');
plot(2*N+1:3*N,y2(2*N+1:end),'*');
axis([1,3*N -0.05 0.05])

% figure
% scatter3(u(1:N,1),u(1:N,2),u(1:N,3),'*');hold on
% scatter3(u(N+1:2*N,1),u(N+1:2*N,2),u(N+1:2*N,3),'o');
% scatter3(u(2*N+1:3*N,1),u(2*N+1:3*N,2),u(2*N+1:3*N,3),'+');