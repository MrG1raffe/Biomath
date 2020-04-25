%% Настройка параметров
a_size = 1000;
u_size = 1000;
a_max = 7 * (7 / 6) ^ 6;
a = linspace(0, a_max, a_size + 1);
a(1) = [];

%% Мультик
nFrames = 500;
mov(1:nFrames) = struct('cdata', [], 'colormap', []);
figure('Name', '(u, f) plot');
for i=1:nFrames
    ind = 2 * i;
    u = linspace(0, a(ind)^(-1/3), u_size);
    f1 = abs(f(u, a(ind)));
    f2 = abs(f(f1, a(ind)));
    f3 = abs(f(f2, a(ind)));
    plot(u, f1, 'b', u, f2, 'r', u, f3, 'g', u, u, 'black--');
    legend('f', 'f^2', 'f^3');
    axis equal;
    %title(['График f(u) при a = ', num2str(a(ind))]);
    %axis([0 1 0 1]);
    mov(i) = getframe;
end
%%
movie(mov)
%%
v = VideoWriter('f_plot2.avi');
open(v);
writeVideo(v, mov);
close(v);
%% График f(u)
ind = ceil(a_size * 0.9 + 1);
u = linspace(0, a(ind)^(-1/3), u_size);
f1 = abs(f(u, a(ind)));
f2 = abs(f(f1, a(ind)));
f3 = f(f2, a(ind));
figure
plot(u, f1, 'b', u, f2, 'r', u, f3, 'g', u, u, 'black--');
legend('f', 'f^2', 'f^3');
%title(['График f(u) при a = ', num2str(a(ind))]);
axis equal;

%% Устойчивость
u_max = a.^(-1/3);
x0 = u_max;
equil_point = fsolve(@(x) stab_fun(x, a), x0);

f_der_eq = zeros(size(a));
for i = 1:numel(a)
    f_der_eq(i) = f_der(equil_point(i), a(i));
end
idx_stab = find(f_der_eq <= -1, 1, 'first');
a_stab = a(idx_stab);
figure('Name', 'Derivative in equilibrium points');
plot(a, f_der_eq, 'b', a, -ones(size(a)), 'black--', a_stab, f_der_eq(idx_stab), 'r*');
disp("Точка потери устойчивости:");
disp(a_stab);

%% Диаграмма равновесий
figure('Name', 'Equilibriums');
plot(a(1:idx_stab-1), equil_point(1:idx_stab-1), 'g', a(idx_stab), equil_point(idx_stab), 'bo', a(idx_stab+1:a_size), equil_point(idx_stab+1:a_size), 'r');
legend('Асимптотически устойчивые равновесия', 'Потеря устойчивости', 'Неустойчивые равновесия');

%% Бифуркационная диаграмма
figure('Name', 'Diagram');
n_diag = 100;
n_stab = 100;
A = zeros(a_size, n_diag);
for i = 1:a_size
     u_cur = 0.1 * a(i)^(-1/3);
     for j = 1:n_stab
         u_cur = f(u_cur, a(i));
     end
     for j = 1:n_diag
         u_cur =  f(u_cur, a(i));
         A(i, j) = u_cur;
     end
end
plot(a, A, '.b');

%% Бифуркационная диаграмма (r = 15-16) устойчивый цикл длины 3 при a = 15.3
figure('Name', 'Diagram');
n_diag = 100;
n_stab = 100;
a_sm = linspace(15, 16, a_size);
A = zeros(a_size, n_diag);
for i = 1:a_size
     u_cur = 0.1 * a_sm(i)^(-1/3);
     for j = 1:n_stab
         u_cur = f(u_cur, a_sm(i));
     end
     for j = 1:n_diag
         u_cur =  f(u_cur, a_sm(i));
         A(i, j) = u_cur;
     end
end
plot(a_sm, A, '.b');

%% Показатель Ляпунова
a_size = 1000;
a_max = 7 * (7 / 6) ^ 6;
a = linspace(0, a_max, a_size + 1);
a(1) = [];
L = zeros(a_size, 1);
for i=1:a_size
    L(i) = lyap(0.1 * a(i)^(-1/3), a(i));
end
plot(a, L, 'b', a, zeros(size(a)), 'black');
%title('Показатель Ляпунова');

%% Цикл длины 3
a = 15.3;
n = 50;
u = zeros(n, 1);
u(1) = 0.3 * a^(-1/3);
for i = 2:n
    u(i) = f(u(i-1), a);
end
plot(u, 's--b', 'MarkerEdgeColor','g');
u1 = f(u(end), a);
u2 = f(u1, a);
u3 = f(u2, a);
disp('Точки цикла длины 3');
disp(u1);
disp(u2);
disp(u3);

%% Пример (устойчивый цикл длины 2)
n = 40;
u = zeros(n, 1);
u(1) = 0.2;
a_cur = 6;
for i=2:n
    u(i) = f(u(i-1), a_cur);
end

plot(u, 's--b', 'MarkerEdgeColor','g');
u1 = f(u(end), a_cur);
u2 = f(u1, a_cur);
disp('Точки цикла длины 2:');
disp(u1);
disp(u2);
%% Пример (все очень хаотично)
n = 100;
u = zeros(n, 1);
u(1) = 0.2;
a_cur = 17;
for i=2:n
    u(i) = f(u(i-1), a_cur);
end

plot(u, 's--b', 'MarkerEdgeColor','g');

%%
n = 200;
u = zeros(n, 1);
u(1) = 0.2;
a_cur = 4;
for i=2:n
    u(i) = f(u(i-1), a_cur);
end

plot(u, 's--b', 'MarkerEdgeColor','g');
%% Задача 2
u_size = 30;
a = 1;
u = linspace(0, a, u_size);
v = u;
[U, V] = meshgrid(u, v);
surf(U, V, func2(U, V, a));
%% Устойчивость
a_size = 400;
a = linspace(0, 1, a_size + 1);
a(1) = [];
x0 = a.^(-1/3);
equil_point = fsolve(@(x) stab_fun(x, a), x0);

d = @(x, a) 0.25 - 12*a*x^2.5;
D = zeros(a_size, 1);
lambda1 = D;
lambda2 = D;
for i=1:a_size
    D(i) = d(equil_point(i), a(i));
    lambda1(i) = 0.25 - 0.5 * sqrt(D(i));
    lambda2(i) = 0.25 + 0.5 * sqrt(D(i));
end
subplot(1, 2, 1);
plot(a, abs(lambda1), 'b');
legend('|\lambda_1|');
subplot(1, 2, 2);
plot(a, abs(lambda2), 'b');
legend('|\lambda_2|');
%% Бифуркационная диаграмма
figure('Name', 'Diagram');
a_size = 1000;
a = linspace(0, 1, a_size);
n_diag = 100;
n_stab = 100;
A = zeros(a_size, n_diag);
for i = 1:a_size
     u_prev = 0.2 * a(i)^(-1/3);
     u_cur = f(u_prev, a(i));
     for j = 1:n_stab
         tmp = u_cur;
         u_cur = func2(u_cur, u_prev, a(i));
         u_prev = tmp;
     end
     for j = 1:n_diag
         tmp = u_cur;
         u_cur = func2(u_cur, u_prev, a(i));
         u_prev = tmp;         
         A(i, j) = u_cur;
     end
end
plot(a, A, '.b');

%% Пример (устойчивость)
a = 1;
n = 80;
u = zeros(n, 1);
u(1) = 0.6* a^(-1/3);
u(2) = f(u(1), a);
for i = 3:n
    u(i) = func2(u(i-1), u(i-2), a);
end
plot(u, 's--b', 'MarkerEdgeColor','g');
