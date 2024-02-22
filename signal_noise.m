%% 计算气溶胶和分子的消光系数
%计算大气分子及气溶胶后向散射系数和消光系数
clc;clear all;close all;
%%
Sa = 39.5 ;           % 气溶胶消光后向散射比  39.5 
Sm = 8*pi/3;          % 大气分子消光后向散射比
delt_z = 15;          % 距离分辨率，m  (15m)
R = 0:delt_z:30000;   % 探测距离区间，m
lambda = 532;         % 发射激光波长，nm
beta_m = 1.54*10^(-6)*(532/lambda)^4.*exp(-R/7000); % 大气分子后向散射系数。(单位:m时)
alpha_m = Sm*beta_m;   %大气分子消光系数
beta_a = (2.47*10^(-6).*exp(-R/2000)+5.13*10^(-9).*exp(-(R-20000).^2/36000000))*(532/lambda); % 气溶胶后向散射系数(单位:m时)
% 原始的公式为下面两个
% beta_m = 1.54*10^(-3)*(532/lambda)^4*exp(-R/7); % 大气分子后向散射系数。(单位:km时？)
% beta_a = (2.47*10^(-3)*exp(-R/2)+5.13*10^(-6)*exp(-(R-20).^2/36))*(532/lambda); % 气溶胶后向散射系数(单位:km时？)
alpha_a = Sa*beta_a;   %气溶胶消光系数
alpha_total_US = alpha_a + alpha_m; % 总的消光系数=气溶胶消光系数加上分子消光系数

intalpha_US = zeros(1,length(alpha_total_US));
for count = 1:length(alpha_total_US)
    int = 0;
    for count1 = 1:count
        int = int + alpha_total_US(count1)*delt_z;
    end
    intalpha_US(count) = int;
end
%% 画图
%%%%%%%%%%%%%%%%分子和气溶胶的后向散射系数图%%%%%%%%%%%%%%%%
figure
semilogx(beta_m,R*1e-3,'b-','linewidth',2,'markersize',25)
hold on
semilogx(beta_a,R*1e-3,'r-','linewidth',2,'markersize',25)
xlabel('Backscatter Coefficient/m^{-1}sr^{-1}','Fontsize',15,'fontname','Times','fontweight','bold');
ylabel('Altitude /km','Fontsize',15,'fontname','Times','fontweight','bold');
set(gca,'FontName','Times New Roman','FontSize',12,'fontweight','bold')
grid on
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);% 这句会出错
legend('molecular','aerosol');

%%%%%%%%%%%%%%%%分子和气溶胶的消光系数图%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%大气总的（分子+气溶胶）消光系数图%%%%%%%%%%%%%%%%%%%%%%
figure
semilogx(alpha_total_US,R*1e-3,'r-','linewidth',2,'markersize',25)
xlabel('Extinciton Coefficient/m^{-1}','Fontsize',15,'fontname','Times','fontweight','bold');
ylabel('Altitude /km','Fontsize',15,'fontname','Times','fontweight','bold');
set(gca,'FontName','Times New Roman','FontSize',12,'fontweight','bold')
grid on
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);

%%%%%%%%%%%%%%%%大气总的积分消光系数图%%%%%%%%%%%%%%%%%%%%%%
figure
plot(intalpha_US,R*1e-3,'r','linewidth',2,'markersize',25)
xlabel('Extinction Coefficient Integral','Fontsize',15,'fontname','Times','fontweight','bold');
ylabel('Altitude /km','Fontsize',15,'fontname','Times','fontweight','bold');
set(gca,'FontName','Times New Roman','FontSize',12,'fontweight','bold')
grid on
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);

%% 单脉冲回波光子
c = 3e8;                 % 光速m/s
t = 10e-9;               % 激光脉宽
fib_ita = 0.2;           % 光纤耦合效率
fil_ita = 0.36*0.8*0.8*0.8*0.8*0.76;          % 整体光学透过率（望远镜+光纤准直+滤光片+扩束+反射镜）
fp_ita = 0.8;            % FP峰值透过率
op_ita = fib_ita*fil_ita*fp_ita; % 光学效率
pc_ita = 0.2;            % PMT量子效率
oe = 1/16;               % 面元有效率

ita = op_ita * pc_ita * oe;   % 系统效率
S_energy = 0.05;         % 单脉冲激光能量 J
Pt = S_energy/t;         % 激光峰值功率 W
r = 0.40;                % 接收望远镜半径m
A = pi*r^2;              % 望远镜接收面积*Oc_rate
h = 6.6262*10^(-34);     % 普朗克常量
freq = 100;              % 重复频率 hz

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
%% 画图    
figure;
semilogy(R*1e-3,Ps ,'b-','linewidth',2,'markersize',25); %单脉冲光子数
ylabel('Photon Count','Fontsize',15,'fontname','Times','fontweight','bold');
xlabel('Altitude /km','Fontsize',15,'fontname','Times','fontweight','bold');
set(gca,'FontName','Times New Roman','FontSize',12,'fontweight','bold')
grid on

%% 时间积分1分钟
                                                             %时间积分回波光子数
% Pn_1min = Pn_ode*100*freq; % 100s时间积分，重复频率为100Hz
int_t=100;%积分时间
M=int_t*freq;% 脉冲个数
Ps_int = Ps*M; % 1min时间积分，重复频率为100Hz
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
% save C:\Users\dell\Desktop\新建文件夹(4)\雪飞的数据\所有高度上都一遍\求我们的系统信噪比\R_US.mat R
% save C:\Users\dell\Desktop\新建文件夹(4)\雪飞的数据\所有高度上都一遍\求我们的系统信噪比\intalpha_US.mat intalpha_US
% save C:\Users\dell\Desktop\新建文件夹(4)\雪飞的数据\所有高度上都一遍\求我们的系统信噪比\photons_single_US.mat photons_single_US
% save G:\瑞利布里渊\大气回波噪声\用背景光求噪声比例因子\求我们的系统信噪比\photons_int.mat photons_int

%% 计算噪声和信噪比
sample_t=100e-9;         % 采集卡的采集时间（时间分辨率）
omega=1e-3;              % 望远镜接收视场角
lambda=355e-9;           % 波长
E=0.26e9;                % 辐射亮度
B_L=1e-9;                % 滤光片带宽
N_b=pc_ita*lambda*E*pi*omega^2*B_L*A*op_ita*sample_t/(h*c*4); %背景光回波光子数
cps=3000;
N_d=cps*sample_t;        % PMT暗计数
noise_single_day=Ps+N_b+N_d; % 白天单脉冲噪声光子数
noise_single_night=Ps+N_d; % 夜晚单脉冲噪声光子数
% noise_int_day=Ps_int+N_b+N_d;% 白天积分脉冲回波光子数
% noise_int_night=Ps_int+N_d;    % 夜晚积分脉冲回波光子数
NSF=1.3598;              % 噪声比例因子

%% 求噪声比例因子
% miu=N_b; %期望
% sigma=20;%标准差
% a=1;
% b=500;
% R = random('Normal',miu,sigma,a,b); %生成期望为0,标准差为 10的(1行500列)500个正态随机数
% NSF=std(R)/sqrt(mean(R));
%% 信噪比
SNR_single_day=Ps./(NSF*sqrt(noise_single_day));
SNR_single_night=Ps./(NSF*sqrt(noise_single_night));
SNR_int_day=Ps*sqrt(M)./(NSF*sqrt(noise_single_day));
SNR_int_night=Ps*sqrt(M)./(NSF*sqrt(noise_single_night));

%% 有效探测高度SNR=1,SNR=10
single_day_1_10=interp1(SNR_single_day(2:end),R(2:end),[1,10],'spline');
single_night_1_10=interp1(SNR_single_night(2:end),R(2:end),[1,10],'spline');
int_day_1_10=interp1(SNR_int_day(2:end),R(2:end),[1,10],'spline');
int_night_1_10=interp1(SNR_int_night(2:end),R(2:end),[1,10],'spline');

%% 画图
figure 
semilogy(R*1e-3,SNR_single_day,'b-','linewidth',2.5,'markersize',25); % 单脉冲白天夜晚信噪比
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
semilogy(R*1e-3,SNR_int_day,'m-','linewidth',2.5,'markersize',25);% 脉冲积分白天夜晚信噪比
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
semilogy(R*1e-3,SNR_single_night,'r','linewidth',2.5,'markersize',25); % 夜晚的单脉冲和积分脉冲信噪比
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
semilogy(R*1e-3,SNR_single_day,'b','linewidth',2.5,'markersize',25); % 白天的单脉冲和积分脉冲信噪比
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
%% 另外一种定义信噪比的公式
SNR1=Ps./sqrt(Ps+2*(N_b+N_d));
SNR2=Ps./sqrt(Ps+2*(0+N_d));
SNR3=Ps*sqrt(M)./sqrt(Ps+2*(N_b+N_d));
SNR4=Ps*sqrt(M)./sqrt(Ps+2*(0+N_d));
SNR1-SNR_single_day;

figure 
semilogy(R*1e-3,SNR_single_day,'b-','linewidth',2.5,'markersize',25); % 单脉冲白天信噪比
hold on
semilogy(R*1e-3,SNR1,'r--','linewidth',2.5,'markersize',25);
ylabel('SNR','Fontsize',15,'fontname','Times','fontweight','bold');
xlabel('Altitude /km','Fontsize',15,'fontname','Times','fontweight','bold'); 
set(gca,'FontName','Times New Roman','FontSize',12,'fontweight','bold')
grid on
legend('method1','method2')

figure 
semilogy(R*1e-3,SNR_single_night,'b-','linewidth',2.5,'markersize',25); % 单脉冲夜晚信噪比
hold on
semilogy(R*1e-3,SNR2,'r--','linewidth',2.5,'markersize',25);
ylabel('SNR','Fontsize',15,'fontname','Times','fontweight','bold');
xlabel('Altitude /km','Fontsize',15,'fontname','Times','fontweight','bold'); 
set(gca,'FontName','Times New Roman','FontSize',12,'fontweight','bold')
grid on
legend('method1','method2')

figure 
semilogy(R*1e-3,SNR_int_day,'b-','linewidth',2.5,'markersize',25); % 1万单脉冲白天信噪比
hold on
semilogy(R*1e-3,SNR3,'r--','linewidth',2.5,'markersize',25);
ylabel('SNR','Fontsize',15,'fontname','Times','fontweight','bold');
xlabel('Altitude /km','Fontsize',15,'fontname','Times','fontweight','bold'); 
set(gca,'FontName','Times New Roman','FontSize',12,'fontweight','bold')
grid on
legend('method1','method2')

figure 
semilogy(R*1e-3,SNR_int_night,'b-','linewidth',2.5,'markersize',25); % 1万脉冲夜晚信噪比
hold on
semilogy(R*1e-3,SNR4,'r--','linewidth',2.5,'markersize',25);
ylabel('SNR','Fontsize',15,'fontname','Times','fontweight','bold');
xlabel('Altitude /km','Fontsize',15,'fontname','Times','fontweight','bold'); 
set(gca,'FontName','Times New Roman','FontSize',12,'fontweight','bold')
grid on
legend('method1','method2')


