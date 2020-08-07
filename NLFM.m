clc,clear,close all
%% Initialization parameters
t_nlfm = 1.3e-6;          %NLFM信号脉宽
b_nlfm = 20e+6;           %NLFM信号带宽
ku = b_nlfm/t_nlfm;       %噪声调制系数
t = t_nlfm;               %被调制信号脉宽
b = 25e+6;                %被调制信号带宽
k = b/t;                  %被调制信号的调频斜率
fs = 400e+6;              %采样频率
fc = 100e+6;              %载频

rng(10);                  %每次运行程序，生成的随机数一定
alp = 2*pi*rand(1)-pi;    %相位随机变量

delt = 1/fs:1/fs:t_nlfm;  %采样点
rng(5);
u = normrnd(0,1,1,length(delt));   %高斯白噪声

%% Build transmitted signal
s_nlfm = zeros(size(delt));
%噪声线性调频信号时域表达式
for i = 1:length(delt)
    s_nlfm(i) = 1/sqrt(t_nlfm)*exp(1i*pi*k*delt(i)^2)*...
        exp(1i*alp)*exp(1i*2*pi*ku*sum(u(1:i)));
end

freq = linspace(0,b_nlfm,length(delt))/1e+6;
s_nlfm_fft = fftshift(fft(s_nlfm));  %频域表达式

%save('NLFM.mat','s_nlfm','s_nlfm_fft');
%
figure(1)
plot(delt,s_nlfm/max(s_nlfm));
xlabel('时间/s'),ylabel('归一化幅度值')
title('NLFM信号时域图')
axis tight

figure(2)
plot(freq,abs(s_nlfm_fft));
xlabel('频率/MHz'),ylabel('幅度值')
title('NLFM信号频域图')
axis tight
