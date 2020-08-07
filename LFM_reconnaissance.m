clc,clear,close all
%%
C = 3e+8;                   %����
T = 10e-6;                  %ʱ��
B = 30e+6;                  %����
K = B/T;                    %��Ƶб��
Fs = 5*B;                   %����Ƶ��
Ts = 1/Fs;                  %�������
R = [9000 10000 10200];     %Ŀ��λ��
RCS = [1 1 1];              %�״�ɢ�����
Rmin = 8500; Rmax = 11500;  
Rwid = Rmax - Rmin;         %ɨ�跶Χ
Twid = 2*Rwid/C;            %��Ų�������ʱ����
Nwid = ceil(Twid/Ts);       %��������

%%
t = linspace(2*Rmin/C,2*Rmax/C,Nwid);     %������
M = length(R);
td = ones(M,1)*t - 2*R'/C * ones(1,Nwid); %???
SNR = [1 0.1 0.01 0.001 10 100 1000];     %���������,�������α�ǿ
for i = 1:7
    Srt1 = RCS*(exp(1i*pi*K*td.^2).*(abs(td)<T/2));                 %�����������״﷢�����Ե�Ƶ�ź�
    n = sqrt(0.5*SNR(i))*(randn(size(Srt1))+1i*randn(size(Srt1)));  %�����ź�
    Srt = Srt1+n;                        %�����д������ź�
    Nchirp=ceil(T/Ts);
    Nfft=2^nextpow2(Nwid+Nwid-1);        %Ϊ�˼ӿ�fft�����ٶ�
    Srw=fft(Srt, Nfft);                  %���������״��ź�Ƶ����
    Srw1=fft(Srt1,Nfft);                 %�����������״��ź�Ƶ����
    t0=linspace(-T/2,T/2,Nchirp);
    St=exp(1i*pi*K*t0.^2);               %�����˾����Ե�Ƶ�ź�
    Sw=fft(St,Nfft);                     %LFMƵ����
    Sot=fftshift(ifft(Srw.*conj(Sw)));   %ʱ�����൱��Ƶ��˻�
                                         %�ڸ���Ҷ���任��ʱ��
    Sot1=fftshift(ifft(Srw1.*conj(Sw))); 
    N0=Nfft/2-Nchirp/2;
    Z=abs(Sot(N0:N0+Nwid-1));
    SNR_dB = -1*10*log10(SNR(i));
    figure
    subplot(211)
    plot(t*1e+6,real(Srt));
    axis tight
    xlabel('us' );ylabel('����');
    title(['�������Ե�Ƶ�ź�ѹ��ǰ,SNR =',num2str(SNR_dB)]);
    subplot(212)
    plot(t*C/2,Z)
    xlabel('Range in meter');ylabel('����');
    title(['�������Ե�Ƶ�ź�ѹ����,SNR =',num2str(SNR_dB)]);
end
