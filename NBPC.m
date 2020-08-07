clc,clear,close all

%% Initialization parameters
tau = 0.1e-6;             %单个码元脉宽
Barker = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];%十三位巴克码
t_nbpc = length(Barker)*tau;             %NBPC信号脉宽
b_nbpc = 10e+6;           %NBPC信号带宽
ku = b_nbpc/t_nbpc;       %噪声调制系数
t = t_nbpc;               %被调制信号脉宽
b = 14e+6;                %被调制信号带宽
fs = 400e+6;              %采样频率
fc = 100e+6;              %载频

rng(10);                  %每次运行程序，生成的随机数一定
alp = 2*pi*rand(1)-pi;    %相位随机变量

delt = 1/fs:1/fs:tau;     %采样点
N = length(delt);         %单个码长采样点的个数
rng(5);
u = normrnd(0,1,1,N*length(Barker));%高斯白噪声

%% Build transmitted signal
s_nbpc0 = zeros(1,length(Barker)*N);

for i=1:length(Barker)
     if Barker(i)==1
        pha = 1;
     else
        pha = -1;
     end
     s_nbpc0(1,(i-1)*N+1:i*N) = pha*1/sqrt(length(Barker)*tau);
end 

delt1 = 1/fs:1/fs:t_nbpc;
s_nbpc = zeros(size(s_nbpc0));
for j = 1:length(delt1)
    s_nbpc(j) = s_nbpc0(j)*exp(1i*2*pi*ku*sum(u(1:j)))*exp(1i*alp);
end

freq = linspace(0,b_nbpc,length(delt1))/1e+6;
s_nbpc_fft = fftshift(fft(s_nbpc)); 

%save('NBPC.mat','s_nbpc','s_nbpc_fft');
% 
% figure(1)
% plot(delt1,s_nbpc/max(s_nbpc));
% xlabel('时间/s'),ylabel('归一化幅度值')
% title('NBPC信号时域图')
% axis tight
% 
% figure(2)
% plot(freq,abs(s_nbpc_fft));
% xlabel('频率/MHz'),ylabel('幅度值')
% title('NBPC信号频域图')
% axis tight
