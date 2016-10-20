%построение графика s(t)=Acos(2pi*f0*t-theta) 0<t<T
clear all %clear workspace
clc %claer command window
close all %закрытие графич.окон
nfig=1; %номер графического окна
%задание основных параметров
T = 4e-3; %0.004
A = 2;
f0 = 1500; %Hz
theta = 3*pi/8;
Ns = 12;%число отсчетов на длительности периода несущей частоты
dt = (1/f0)/Ns;%шаг по оси t
t = 0:dt:T;
%вычисление
s = A*cos(2*pi*f0*t-theta);
%построение
figure(nfig); nfig = nfig+1;
plot(t,s,'-bo', 'LineWidth',2);
%   <эксперименты>
hold on
Ns = 32;
N = 20;
dt = (1/f0)/Ns;%шаг по оси t
t = 0:dt:T;
dt1 = (1/f0)/N;%шаг по оси t
t1 = 0:dt1:T;
%вычисление
s = A*cos(2*pi*f0*t-theta-pi);
s1 = A*cos(2*pi*f0*t1-theta-pi/2);
plot(t,s,'-mp', t1, s1, '-r*');
hold off
%   </эксперименты>

%эти параметры(color/linetype argument) в любом порядке:
%буква o ставит кружочки на места точек
% "p"/pentagram звездочки
% "h"/hexagram звездочки
% "d"/diamond ставит ромбики на места точек
%"+" плюсы на место точек
% "-" соединяет точки в линию
% "-."/":"/"--" пунктирная линия
%больше - в Ctrl + D
%цвета blue, red, yellow, green, magenta, white, cyan, blacK
grid on
xlabel('t')
ylabel('s(t)')
axis([0-10*dt, T+10*dt, -1.5*A, 1.5*A])