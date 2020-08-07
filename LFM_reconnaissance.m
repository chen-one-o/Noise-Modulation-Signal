clc,clear,close all
%%
C = 3e+8;                   %光速
T = 10e-6;                  %时宽
B = 30e+6;                  %带宽
K = B/T;                    %调频斜率
Fs = 5*B;                   %抽样频率
Ts = 1/Fs;                  %抽样间隔
R = [9000 10000 10200];     %目标位置
RCS = [1 1 1];              %雷达散射截面
Rmin = 8500; Rmax = 11500;  
Rwid = Rmax - Rmin;         %扫描范围
Twid = 2*Rwid/C;            %电磁波往返的时间间隔
Nwid = ceil(Twid/Ts);       %抽样点数

%%
t = linspace(2*Rmin/C,2*Rmax/C,Nwid);     %采样点
M = length(R);
td = ones(M,1)*t - 2*R'/C * ones(1,Nwid); %???
SNR = [1 0.1 0.01 0.001 10 100 1000];     %信噪比向量,噪声依次变强
for i = 1:7
    Srt1 = RCS*(exp(1i*pi*K*td.^2).*(abs(td)<T/2));                 %满足条件的雷达发射线性调频信号
    n = sqrt(0.5*SNR(i))*(randn(size(Srt1))+1i*randn(size(Srt1)));  %噪声信号
    Srt = Srt1+n;                        %大气中传播的信号
    Nchirp=ceil(T/Ts);
    Nfft=2^nextpow2(Nwid+Nwid-1);        %为了加快fft计算速度
    Srw=fft(Srt, Nfft);                  %加噪声的雷达信号频域表达
    Srw1=fft(Srt1,Nfft);                 %不加噪声的雷达信号频域表达
    t0=linspace(-T/2,T/2,Nchirp);
    St=exp(1i*pi*K*t0.^2);               %正儿八经线性调频信号
    Sw=fft(St,Nfft);                     %LFM频域表达
    Sot=fftshift(ifft(Srw.*conj(Sw)));   %时域卷积相当于频域乘积
                                         %在傅里叶反变换到时域
    Sot1=fftshift(ifft(Srw1.*conj(Sw))); 
    N0=Nfft/2-Nchirp/2;
    Z=abs(Sot(N0:N0+Nwid-1));
    SNR_dB = -1*10*log10(SNR(i));
    figure
    subplot(211)
    plot(t*1e+6,real(Srt));
    axis tight
    xlabel('us' );ylabel('幅度');
    title(['加噪线性调频信号压缩前,SNR =',num2str(SNR_dB)]);
    subplot(212)
    plot(t*C/2,Z)
    xlabel('Range in meter');ylabel('幅度');
    title(['加噪线性调频信号压缩后,SNR =',num2str(SNR_dB)]);
end
