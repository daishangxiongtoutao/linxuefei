%% �������ܽ��ͷ��ӵ�����ϵ��
%����������Ӽ����ܽ�����ɢ��ϵ��������ϵ��
clc;clear all;close all;
%%
Sa = 39.5 ;           % ���ܽ��������ɢ���  39.5 
Sm = 8*pi/3;          % ���������������ɢ���
delt_z = 15;          % ����ֱ��ʣ�m  (15m)
R = 0:delt_z:30000;   % ̽��������䣬m
lambda = 532;         % ���伤�Ⲩ����nm
beta_m = 1.54*10^(-6)*(532/lambda)^4.*exp(-R/7000); % �������Ӻ���ɢ��ϵ����(��λ:mʱ)
alpha_m = Sm*beta_m;   %������������ϵ��
beta_a = (2.47*10^(-6).*exp(-R/2000)+5.13*10^(-9).*exp(-(R-20000).^2/36000000))*(532/lambda); % ���ܽ�����ɢ��ϵ��(��λ:mʱ)
% ԭʼ�Ĺ�ʽΪ��������
% beta_m = 1.54*10^(-3)*(532/lambda)^4*exp(-R/7); % �������Ӻ���ɢ��ϵ����(��λ:kmʱ��)
% beta_a = (2.47*10^(-3)*exp(-R/2)+5.13*10^(-6)*exp(-(R-20).^2/36))*(532/lambda); % ���ܽ�����ɢ��ϵ��(��λ:kmʱ��)
alpha_a = Sa*beta_a;   %���ܽ�����ϵ��
alpha_total_US = alpha_a + alpha_m; % �ܵ�����ϵ��=���ܽ�����ϵ�����Ϸ�������ϵ��

intalpha_US = zeros(1,length(alpha_total_US));
for count = 1:length(alpha_total_US)
    int = 0;
    for count1 = 1:count
        int = int + alpha_total_US(count1)*delt_z;
    end
    intalpha_US(count) = int;
end
%% ��ͼ
%%%%%%%%%%%%%%%%���Ӻ����ܽ��ĺ���ɢ��ϵ��ͼ%%%%%%%%%%%%%%%%
figure
semilogx(beta_m,R*1e-3,'b-','linewidth',2,'markersize',25)
hold on
semilogx(beta_a,R*1e-3,'r-','linewidth',2,'markersize',25)
xlabel('Backscatter Coefficient/m^{-1}sr^{-1}','Fontsize',15,'fontname','Times','fontweight','bold');
ylabel('Altitude /km','Fontsize',15,'fontname','Times','fontweight','bold');
set(gca,'FontName','Times New Roman','FontSize',12,'fontweight','bold')
grid on
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);% �������
legend('molecular','aerosol');

%%%%%%%%%%%%%%%%���Ӻ����ܽ�������ϵ��ͼ%%%%%%%%%%%%%%%%%%%%%%
figure
semilogx(alpha_m,R*1e-3,'b-','linewidth',2,'markersize',25)
hold on
semilogx(alpha_a,R*1e-3,'r-','linewidth',2,'markersize',25)
xlabel('Extinciton Coefficient/m^{-1}','Fontsize',15,'fontname','Times','fontweight','bold');
ylabel('Altitude /km','Fontsize',15,'fontname','Times','fontweight','bold');
set(gca,'FontName','Times New Roman','FontSize',12,'fontweight','bold')
grid on
% set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
legend('molecular','aerosol');

%%%%%%%%%%%%%%%%�����ܵģ�����+���ܽ�������ϵ��ͼ%%%%%%%%%%%%%%%%%%%%%%
figure
semilogx(alpha_total_US,R*1e-3,'r-','linewidth',2,'markersize',25)
xlabel('Extinciton Coefficient/m^{-1}','Fontsize',15,'fontname','Times','fontweight','bold');
ylabel('Altitude /km','Fontsize',15,'fontname','Times','fontweight','bold');
set(gca,'FontName','Times New Roman','FontSize',12,'fontweight','bold')
grid on
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);

%%%%%%%%%%%%%%%%�����ܵĻ�������ϵ��ͼ%%%%%%%%%%%%%%%%%%%%%%
figure
plot(intalpha_US,R*1e-3,'r','linewidth',2,'markersize',25)
xlabel('Extinction Coefficient Integral','Fontsize',15,'fontname','Times','fontweight','bold');
ylabel('Altitude /km','Fontsize',15,'fontname','Times','fontweight','bold');
set(gca,'FontName','Times New Roman','FontSize',12,'fontweight','bold')
grid on
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);

%% ������ز�����
c = 3e8;                 % ����m/s
t = 10e-9;               % ��������
fib_ita = 0.2;           % �������Ч��
fil_ita = 0.36*0.8*0.8*0.8*0.8*0.76;          % �����ѧ͸���ʣ���Զ��+����׼ֱ+�˹�Ƭ+����+���侵��
fp_ita = 0.8;            % FP��ֵ͸����
op_ita = fib_ita*fil_ita*fp_ita; % ��ѧЧ��
pc_ita = 0.2;            % PMT����Ч��
oe = 1/16;               % ��Ԫ��Ч��

ita = op_ita * pc_ita * oe;   % ϵͳЧ��
S_energy = 0.05;         % �����弤������ J
Pt = S_energy/t;         % �����ֵ���� W
r = 0.40;                % ������Զ���뾶m
A = pi*r^2;              % ��Զ���������*Oc_rate
h = 6.6262*10^(-34);     % ���ʿ˳���
freq = 100;              % �ظ�Ƶ�� hz

sys_const = (S_energy*lambda*10^-9)*A* ita*delt_z/(h*c);
for i =1:length(R)
    Pn_a1(i) = sys_const*beta_a(i)/(R(i)^2);
    Pn_a2(i) = exp(-2*intalpha_US(i));
    Pn_a(i) = Pn_a1(i)*Pn_a2(i);
    Pn_m1(i) = sys_const*beta_m(i)/(R(i)^2);
    Pn_m2(i) = exp(-2*intalpha_US(i));
    Pn_m(i) = Pn_m1(i)*Pn_m2(i);
end
    Ps = (Pn_m + Pn_a);                                        
%% ��ͼ    
figure;
semilogy(R*1e-3,Ps ,'b-','linewidth',2,'markersize',25); %�����������
ylabel('Photon Count','Fontsize',15,'fontname','Times','fontweight','bold');
xlabel('Altitude /km','Fontsize',15,'fontname','Times','fontweight','bold');
set(gca,'FontName','Times New Roman','FontSize',12,'fontweight','bold')
grid on

%% ʱ�����1����
                                                             %ʱ����ֻز�������
% Pn_1min = Pn_ode*100*freq; % 100sʱ����֣��ظ�Ƶ��Ϊ100Hz
int_t=100;%����ʱ��
M=int_t*freq;% �������
Ps_int = Ps*M; % 1minʱ����֣��ظ�Ƶ��Ϊ100Hz
figure 
semilogy(R*1e-3,Ps_int,'m','linewidth',1.5,'markersize',25);
ylabel('Photon Count','Fontsize',15,'fontname','Times','fontweight','bold');
xlabel('Altitude /km','Fontsize',15,'fontname','Times','fontweight','bold');
set(gca,'FontName','Times New Roman','FontSize',12,'fontweight','bold')
grid on
hold off

%%
% photons_single_US=Ps ;
% photons_int=Ps_int;
% save C:\Users\dell\Desktop\�½��ļ���(4)\ѩ�ɵ�����\���и߶��϶�һ��\�����ǵ�ϵͳ�����\R_US.mat R
% save C:\Users\dell\Desktop\�½��ļ���(4)\ѩ�ɵ�����\���и߶��϶�һ��\�����ǵ�ϵͳ�����\intalpha_US.mat intalpha_US
% save C:\Users\dell\Desktop\�½��ļ���(4)\ѩ�ɵ�����\���и߶��϶�һ��\�����ǵ�ϵͳ�����\photons_single_US.mat photons_single_US
% save G:\��������Ԩ\�����ز�����\�ñ�������������������\�����ǵ�ϵͳ�����\photons_int.mat photons_int

%% ���������������
sample_t=100e-9;         % �ɼ����Ĳɼ�ʱ�䣨ʱ��ֱ��ʣ�
omega=1e-3;              % ��Զ�������ӳ���
lambda=355e-9;           % ����
E=0.26e9;                % ��������
B_L=1e-9;                % �˹�Ƭ����
N_b=pc_ita*lambda*E*pi*omega^2*B_L*A*op_ita*sample_t/(h*c*4); %������ز�������
cps=3000;
N_d=cps*sample_t;        % PMT������
noise_single_day=Ps+N_b+N_d; % ���쵥��������������
noise_single_night=Ps+N_d; % ҹ����������������
% noise_int_day=Ps_int+N_b+N_d;% �����������ز�������
% noise_int_night=Ps_int+N_d;    % ҹ���������ز�������
NSF=1.3598;              % ������������

%% ��������������
% miu=N_b; %����
% sigma=20;%��׼��
% a=1;
% b=500;
% R = random('Normal',miu,sigma,a,b); %��������Ϊ0,��׼��Ϊ 10��(1��500��)500����̬�����
% NSF=std(R)/sqrt(mean(R));
%% �����
SNR_single_day=Ps./(NSF*sqrt(noise_single_day));
SNR_single_night=Ps./(NSF*sqrt(noise_single_night));
SNR_int_day=Ps*sqrt(M)./(NSF*sqrt(noise_single_day));
SNR_int_night=Ps*sqrt(M)./(NSF*sqrt(noise_single_night));

%% ��Ч̽��߶�SNR=1,SNR=10
single_day_1_10=interp1(SNR_single_day(2:end),R(2:end),[1,10],'spline');
single_night_1_10=interp1(SNR_single_night(2:end),R(2:end),[1,10],'spline');
int_day_1_10=interp1(SNR_int_day(2:end),R(2:end),[1,10],'spline');
int_night_1_10=interp1(SNR_int_night(2:end),R(2:end),[1,10],'spline');

%% ��ͼ
figure 
semilogy(R*1e-3,SNR_single_day,'b-','linewidth',2.5,'markersize',25); % ���������ҹ�������
hold on
semilogy(R*1e-3,SNR_single_night,'r-','linewidth',2.5,'markersize',25);
hold on
line([0,50],[1,1],'color','y','linestyle',':','LineWidth',2);
line([0,50],[10,10],'color','c','linestyle',':','LineWidth',2);
ylabel('SNR','Fontsize',15,'fontname','Times','fontweight','bold');
xlabel('Altitude /km','Fontsize',15,'fontname','Times','fontweight','bold'); 
set(gca,'FontName','Times New Roman','FontSize',12,'fontweight','bold')
grid on
legend('day','night','SNR=1','SNR=10')

figure 
semilogy(R*1e-3,SNR_int_day,'m-','linewidth',2.5,'markersize',25);% ������ְ���ҹ�������
hold on
semilogy(R*1e-3,SNR_int_night,'k-','linewidth',2.5,'markersize',25);
hold on
line([0,50],[1,1],'color','y','linestyle',':','LineWidth',2);
line([0,50],[10,10],'color','c','linestyle',':','LineWidth',2);
ylabel('SNR','Fontsize',15,'fontname','Times','fontweight','bold');
xlabel('Altitude /km','Fontsize',15,'fontname','Times','fontweight','bold');
set(gca,'FontName','Times New Roman','FontSize',12,'fontweight','bold')
grid on
legend('day','night','SNR=1','SNR=10')

figure 
semilogy(R*1e-3,SNR_single_night,'r','linewidth',2.5,'markersize',25); % ҹ��ĵ�����ͻ������������
hold on
semilogy(R*1e-3,SNR_int_night,'k','linewidth',2.5,'markersize',25);
hold on
line([0,50],[1,1],'color','y','linestyle',':','LineWidth',2);
line([0,50],[10,10],'color','c','linestyle',':','LineWidth',2);
ylabel('SNR','Fontsize',15,'fontname','Times','fontweight','bold');
xlabel('Altitude /km','Fontsize',15,'fontname','Times','fontweight','bold');
set(gca,'FontName','Times New Roman','FontSize',12,'fontweight','bold')
grid on
legend('pulse=1','pulse=10000','SNR=1','SNR=10')

figure 
semilogy(R*1e-3,SNR_single_day,'b','linewidth',2.5,'markersize',25); % ����ĵ�����ͻ������������
hold on
semilogy(R*1e-3,SNR_int_day,'m','linewidth',2.5,'markersize',25);
hold on
line([0,50],[1,1],'color','y','linestyle',':','LineWidth',2);
line([0,50],[10,10],'color','c','linestyle',':','LineWidth',2);
ylabel('SNR','Fontsize',15,'fontname','Times','fontweight','bold');
xlabel('Altitude /km','Fontsize',15,'fontname','Times','fontweight','bold');
set(gca,'FontName','Times New Roman','FontSize',12,'fontweight','bold')
grid on
legend('pulse=1','pulse=10000','SNR=1','SNR=10')
%% ����һ�ֶ�������ȵĹ�ʽ
SNR1=Ps./sqrt(Ps+2*(N_b+N_d));
SNR2=Ps./sqrt(Ps+2*(0+N_d));
SNR3=Ps*sqrt(M)./sqrt(Ps+2*(N_b+N_d));
SNR4=Ps*sqrt(M)./sqrt(Ps+2*(0+N_d));
SNR1-SNR_single_day;

figure 
semilogy(R*1e-3,SNR_single_day,'b-','linewidth',2.5,'markersize',25); % ��������������
hold on
semilogy(R*1e-3,SNR1,'r--','linewidth',2.5,'markersize',25);
ylabel('SNR','Fontsize',15,'fontname','Times','fontweight','bold');
xlabel('Altitude /km','Fontsize',15,'fontname','Times','fontweight','bold'); 
set(gca,'FontName','Times New Roman','FontSize',12,'fontweight','bold')
grid on
legend('method1','method2')

figure 
semilogy(R*1e-3,SNR_single_night,'b-','linewidth',2.5,'markersize',25); % ������ҹ�������
hold on
semilogy(R*1e-3,SNR2,'r--','linewidth',2.5,'markersize',25);
ylabel('SNR','Fontsize',15,'fontname','Times','fontweight','bold');
xlabel('Altitude /km','Fontsize',15,'fontname','Times','fontweight','bold'); 
set(gca,'FontName','Times New Roman','FontSize',12,'fontweight','bold')
grid on
legend('method1','method2')

figure 
semilogy(R*1e-3,SNR_int_day,'b-','linewidth',2.5,'markersize',25); % 1��������������
hold on
semilogy(R*1e-3,SNR3,'r--','linewidth',2.5,'markersize',25);
ylabel('SNR','Fontsize',15,'fontname','Times','fontweight','bold');
xlabel('Altitude /km','Fontsize',15,'fontname','Times','fontweight','bold'); 
set(gca,'FontName','Times New Roman','FontSize',12,'fontweight','bold')
grid on
legend('method1','method2')

figure 
semilogy(R*1e-3,SNR_int_night,'b-','linewidth',2.5,'markersize',25); % 1������ҹ�������
hold on
semilogy(R*1e-3,SNR4,'r--','linewidth',2.5,'markersize',25);
ylabel('SNR','Fontsize',15,'fontname','Times','fontweight','bold');
xlabel('Altitude /km','Fontsize',15,'fontname','Times','fontweight','bold'); 
set(gca,'FontName','Times New Roman','FontSize',12,'fontweight','bold')
grid on
legend('method1','method2')


