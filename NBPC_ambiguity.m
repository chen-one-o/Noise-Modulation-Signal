clc,clear,close all;

%% NBPC ambiguity figure
load NBPC.mat
N = size(s_nbpc,2);%NBPC���г���
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
%��U(f)
V = U;
delay = linspace(-tau,tau,nfft);
freq_del = 12/tau/100;%-6/tau��6/tau�ֳ�100��
vfft = fft(V,nfft);

j = 0;
for freq = -6/tau:freq_del:6/tau %Ƶ����ȡ101����
    j = j+1;
    exf = exp(sqrt(-1)*2.*pi.*freq.*delay);%e^(j2*pi*fd*delay)
    u_times_exf = U .* exf;%s(t)*e^(j2*pi*fd*delay),ΪƵ����׼��
    ufft = fft(u_times_exf,nfft);%S(f-fd)
    prob = ufft.*conj(vfft);% S(f-fd)*����[S(f)]
    ambig(:,j) = fftshift(abs(ifft(prob))');%ģ������ = |���ֺ� ����[S(f)]*S(f-fd)*e^(-j*2*pi*f*tau)|^2
end

%% NLFMģ������ͼ��ģ��ͼ+����ֱ���ͼ+�ٶȷֱ���ͼ��
freq = -6/tau:freq_del:6/tau;
Freq = 10*(1e+6)*freq/(6/tau);
delay = 1.3*linspace(-N,N,nfft)*(1e-6)/N;

figure(1)
mesh(Freq,delay,ambig./max(max(ambig)));%����max����ÿһ�е����ֵ������max������һ�е����ֵ
title('NBPC�ź�ģ������ͼ');
axis tight
xlabel('frequency/MHz');ylabel('delay/\mus');zlabel('ambiguity function');

figure(2)
%value = 10*N;
plot(delay,ambig(:,51)/max(ambig(:,51)));
xlabel('delay/us');ylabel('normalized ambiguity cut for f=0');
title('NBPC�źž���ģ������ͼ');
grid on
axis tight

figure(3)
plot(Freq,ambig(4096,:)/(max(ambig(4096,:))));
xlabel('freqency/MHz');ylabel('normalized ambiguity cut for \tau=0');
title('NBPC�ź��ٶ�ģ������ͼ');
grid on
axis tight