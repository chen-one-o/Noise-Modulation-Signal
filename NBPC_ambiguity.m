clc,clear,close all;

%% NBPC ambiguity figure
load NBPC.mat
N = size(s_nbpc,2);%NBPC序列长度
tau = N;

samp_num = N * 10;
n = ceil(log(samp_num)/log(2));
nfft = 2^n;

U(1:nfft) = 0;
j = 0;
for index = 1:10:samp_num 
    j = j+1;
    U(index:index+10-1) = s_nbpc(j);
end
%求U(f)
V = U;
delay = linspace(-tau,tau,nfft);
freq_del = 12/tau/100;%-6/tau到6/tau分成100份
vfft = fft(V,nfft);

j = 0;
for freq = -6/tau:freq_del:6/tau %频率轴取101个点
    j = j+1;
    exf = exp(sqrt(-1)*2.*pi.*freq.*delay);%e^(j2*pi*fd*delay)
    u_times_exf = U .* exf;%s(t)*e^(j2*pi*fd*delay),为频移做准备
    ufft = fft(u_times_exf,nfft);%S(f-fd)
    prob = ufft.*conj(vfft);% S(f-fd)*共轭[S(f)]
    ambig(:,j) = fftshift(abs(ifft(prob))');%模糊函数 = |积分号 共轭[S(f)]*S(f-fd)*e^(-j*2*pi*f*tau)|^2
end

%% NLFM模糊函数图（模糊图+距离分辨力图+速度分辨力图）
freq = -6/tau:freq_del:6/tau;
Freq = 10*(1e+6)*freq/(6/tau);
delay = 1.3*linspace(-N,N,nfft)*(1e-6)/N;

figure(1)
mesh(Freq,delay,ambig./max(max(ambig)));%里层的max：求每一列的最大值；外层的max：求那一行的最大值
title('NBPC信号模糊函数图');
axis tight
xlabel('frequency/MHz');ylabel('delay/\mus');zlabel('ambiguity function');

figure(2)
%value = 10*N;
plot(delay,ambig(:,51)/max(ambig(:,51)));
xlabel('delay/us');ylabel('normalized ambiguity cut for f=0');
title('NBPC信号距离模糊函数图');
grid on
axis tight

figure(3)
plot(Freq,ambig(4096,:)/(max(ambig(4096,:))));
xlabel('freqency/MHz');ylabel('normalized ambiguity cut for \tau=0');
title('NBPC信号速度模糊函数图');
grid on
axis tight