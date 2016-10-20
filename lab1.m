%���� 1 ������� � �������
clear all %clear workspace
clc %claer command window
close all %�������� ������.����
Vmod = 1200;
Vinf = 4800;
f0 = 1800;
A = 1;
Ns = 32;
%���������� ����������� ����������
T = 1/Vmod;
q = pow2(Vinf*T);
sq = sqrt(q);
dt = (1/f0)/Ns;
t = 0:dt:T;
phi1 = sqrt(2/T)*cos(2*pi*f0*t);
phi2 = sqrt(2/T)*sin(2*pi*f0*t);
%��������� ������ (q �����, t ��������)
s = zeros(q,length(t)); %��������� ������
%  figure plot(t, phi1, '-og', t, phi2, '-oc'); grid on
for i = 0 : q-1
    i2 = mod(i,sq); % i = sq*i1 + i2
    i1 = (i-i2)/sq;
    si1 = A*(1-(2*i1)/(sq-1));
    si2 = A*(1-(2*i2)/(sq-1));
    s(i+1,:) = si1*phi1 + si2*phi2;%��������� ������
end
% ���������� ��������
nfig = 1;
figure(nfig); nfig = nfig+1;
%hold on xlabel('t') ylabel('s(t)')
for i = 1:q
     h(i)=subplot(sq,sq,i);    
     plot(t,s(i,:),'m-','LineWidth',2);
     title(['s_{',int2str(i),'}','(t)']) %si(t)
     grid on
end
% ���������� �� ��� ���� �� 8 �������� - ��� ��� ������ for i = q/2+1:q
%      h(i)=subplot(sq/2,sq,i-q/2); plot(t,s(i,:),'g-','LineWidth',2);
%      title(i) grid on
% end

% ������������� ���������� ��� ��� ���� ��������
axis(h,[0-3*dt T+3*dt -75 75])

% �������
%������� �������������� �����  phi1 � phi2
f = f0-1/T-300 : 100 : f0+1/T+300;%f0-1/T-300
%������� ��� �������� (q �����, f ��������)
S = zeros(q, length(f));
Phi1 = A/2*T*(sinc(pi*T*(f-f0))+sinc(pi*T*(f+f0)));
Phi2 = A/(2*j)*T*(sinc(pi*T*(f-f0))-sinc(pi*T*(f+f0)));
for i = 0 : q-1
    i2 = mod(i,sq); % i = sq*i1 + i2
    i1 = (i-i2)/sq;
    si1 = A*(1-(2*i1)/(sq-1));
    si2 = A*(1-(2*i2)/(sq-1));
    S(i+1,:) = si1*Phi1 + si2*Phi2; %��������� ������
end
figure(nfig); nfig = nfig+1;  %����������
for i = 1 : q
    g(i) = subplot(sq, sq, i);
    plot(f,abs(S(i,:)),'g-','LineWidth',2);
    title(['|S_{',int2str(i),'}|(f)']) % |Si|(f)
    grid on
end
 axis(g, [f(1) f(end) -0.1e-4 6.3e-4])
 
 % ������ ������������������ ��������
 %������ ��������� ������������������� �������� ����� N(������ ������
 % �����) �� ��� �������� s(t)
 % �������� ���-��, ��������� �� ������ � �������� � ���������� ���
 N1 = 7; 
 N2 = 3; 
 N3 = 10;
 %������������� 1-16
 i1 = [2,1,9,1,8,14,7];  %����� 7
 i2 = [5,8,15];
 i3 = [3,3,6,11,12,1,7,5,10,8]; %����� 10

 % S1 - ������� �� ���������*�
 S1 = zeros(1, length(f));
 S2 = zeros(1, length(f));
 S3 = zeros(1, length(f));
 
% ����������������� ��� ������ %
 s1 = []; 
 s2 = [];
 s3 = [];

 for l = 1 : N1
     s1 = cat(2, s1,  s(i1(l),:) );
     S1(1,:) = S1(1,:) + S(i1(l),:).*exp(-1i*2*pi*f*l*T);
 end
for l = 1:N2 
    s2 = cat(2, s2, s(i2(l),:));
    S2(1,:) = S2(1,:) + S(i2(l),:).*exp(-1i*2*pi*f*l*T);
end
for l = 1:N3 
    s3 = cat(2, s3, s(i3(l),:));
    S3(1,:) = S3(1,:) + S(i3(l),:).*exp(-1i*2*pi*f*l*T);
end
 
  t1 = 0:dt:T*length(i1)+6*dt; % ����� ��� ������ ����-�� �� ������
  t2 = 0:dt:T*N2+(N2-1)*dt;
  t3 = 0:dt:T*N3+(N3-1)*dt;
  figure(nfig); nfig = nfig+1;
  plot(t1,s1(1,:),'-r','LineWidth',1);  %���������� ������-�� 1
  title('s_1(t)');
  grid on
  figure(nfig); nfig = nfig+1;
  plot(t2, s2(1,:), '-b','LineWidth',2);  %���������� ������-�� 2
  title('s_2(t)');
  grid on
  figure(nfig); nfig = nfig+1;
  plot(t3, s3(1,:), '-c','LineWidth',2); %���������� ������-�� 3
  title('s_3(t)');
  grid on
  
  figure(nfig); nfig = nfig+1; %���������� �������� ������-���
  hold on
  plot(f,abs(S1),'-r');    %���������� ������� ������-�� 1
  plot(f,abs(S2),'-b','LineWidth',2);    %���������� ������� ������-�� 2
  plot(f, abs(S3), '-c','LineWidth',2)  %���������� ������� ������-�� 3
  hold off
  title('S_{i}(f)');
  grid on
