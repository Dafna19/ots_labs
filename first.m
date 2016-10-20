%���������� ������� s(t)=Acos(2pi*f0*t-theta) 0<t<T
clear all %clear workspace
clc %claer command window
close all %�������� ������.����
nfig=1; %����� ������������ ����
%������� �������� ����������
T = 4e-3; %0.004
A = 2;
f0 = 1500; %Hz
theta = 3*pi/8;
Ns = 12;%����� �������� �� ������������ ������� ������� �������
dt = (1/f0)/Ns;%��� �� ��� t
t = 0:dt:T;
%����������
s = A*cos(2*pi*f0*t-theta);
%����������
figure(nfig); nfig = nfig+1;
plot(t,s,'-bo', 'LineWidth',2);
%   <������������>
hold on
Ns = 32;
N = 20;
dt = (1/f0)/Ns;%��� �� ��� t
t = 0:dt:T;
dt1 = (1/f0)/N;%��� �� ��� t
t1 = 0:dt1:T;
%����������
s = A*cos(2*pi*f0*t-theta-pi);
s1 = A*cos(2*pi*f0*t1-theta-pi/2);
plot(t,s,'-mp', t1, s1, '-r*');
hold off
%   </������������>

%��� ���������(color/linetype argument) � ����� �������:
%����� o ������ �������� �� ����� �����
% "p"/pentagram ���������
% "h"/hexagram ���������
% "d"/diamond ������ ������� �� ����� �����
%"+" ����� �� ����� �����
% "-" ��������� ����� � �����
% "-."/":"/"--" ���������� �����
%������ - � Ctrl + D
%����� blue, red, yellow, green, magenta, white, cyan, blacK
grid on
xlabel('t')
ylabel('s(t)')
axis([0-10*dt, T+10*dt, -1.5*A, 1.5*A])