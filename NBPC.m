clc,clear,close all

%% Initialization parameters
tau = 0.1e-6;             %������Ԫ����
Barker = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];%ʮ��λ�Ϳ���
t_nbpc = length(Barker)*tau;             %NBPC�ź�����
b_nbpc = 10e+6;           %NBPC�źŴ���
ku = b_nbpc/t_nbpc;       %��������ϵ��
t = t_nbpc;               %�������ź�����
b = 14e+6;                %�������źŴ���
fs = 400e+6;              %����Ƶ��
fc = 100e+6;              %��Ƶ

rng(10);                  %ÿ�����г������ɵ������һ��
alp = 2*pi*rand(1)-pi;    %��λ�������

delt = 1/fs:1/fs:tau;     %������
N = length(delt);         %�����볤������ĸ���
rng(5);
u = normrnd(0,1,1,N*length(Barker));%��˹������

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
% xlabel('ʱ��/s'),ylabel('��һ������ֵ')
% title('NBPC�ź�ʱ��ͼ')
% axis tight
% 
% figure(2)
% plot(freq,abs(s_nbpc_fft));
% xlabel('Ƶ��/MHz'),ylabel('����ֵ')
% title('NBPC�ź�Ƶ��ͼ')
% axis tight
