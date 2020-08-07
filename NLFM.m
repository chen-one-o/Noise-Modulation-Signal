clc,clear,close all
%% Initialization parameters
t_nlfm = 1.3e-6;          %NLFM�ź�����
b_nlfm = 20e+6;           %NLFM�źŴ���
ku = b_nlfm/t_nlfm;       %��������ϵ��
t = t_nlfm;               %�������ź�����
b = 25e+6;                %�������źŴ���
k = b/t;                  %�������źŵĵ�Ƶб��
fs = 400e+6;              %����Ƶ��
fc = 100e+6;              %��Ƶ

rng(10);                  %ÿ�����г������ɵ������һ��
alp = 2*pi*rand(1)-pi;    %��λ�������

delt = 1/fs:1/fs:t_nlfm;  %������
rng(5);
u = normrnd(0,1,1,length(delt));   %��˹������

%% Build transmitted signal
s_nlfm = zeros(size(delt));
%�������Ե�Ƶ�ź�ʱ����ʽ
for i = 1:length(delt)
    s_nlfm(i) = 1/sqrt(t_nlfm)*exp(1i*pi*k*delt(i)^2)*...
        exp(1i*alp)*exp(1i*2*pi*ku*sum(u(1:i)));
end

freq = linspace(0,b_nlfm,length(delt))/1e+6;
s_nlfm_fft = fftshift(fft(s_nlfm));  %Ƶ����ʽ

%save('NLFM.mat','s_nlfm','s_nlfm_fft');
%
figure(1)
plot(delt,s_nlfm/max(s_nlfm));
xlabel('ʱ��/s'),ylabel('��һ������ֵ')
title('NLFM�ź�ʱ��ͼ')
axis tight

figure(2)
plot(freq,abs(s_nlfm_fft));
xlabel('Ƶ��/MHz'),ylabel('����ֵ')
title('NLFM�ź�Ƶ��ͼ')
axis tight
