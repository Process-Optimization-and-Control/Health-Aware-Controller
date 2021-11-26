clear all
clc

%load('openLoopAnalysis_v1')
load('openLoopAnalysis')

%% simulation configuration
nMC = 200; 
simLength = 300; % 300 | 450 | 500 | 550

% Initial condition
[x0,z0,u0] = InitialConditionGasLift_5;

% System parameters
par = ParametersGasLift;

%states to measurement mapping function
H = zeros(16,length(z0));
%pai - annulus pressure, well 1-3
H(1,1) = 1; 
H(2,2) = 1; 
H(3,3) = 1; 
%pwh - well head pressure, well+ 1-3
H(4,4) = 1; 
H(5,5) = 1; 
H(6,6) = 1; 
%wro - wellhead gas production rate, well 1-3
H(7,25) = 1; 
H(8,26) = 1; 
H(9,27) = 1; 
%wrg - wellhead oil production rate, well 1-3
H(10,28) = 1; 
H(11,29) = 1; 
H(12,30) = 1; 
%prh - riser head pressure
H(13,37) = 1;  
%pm - manifold pressure
H(14,39) = 1;
%wto - riser head total oil production rate
H(15,41) = 1;
%wtg - riser head total gas production rate
H(16,42) = 1; 
par.H = H;

% generating sand profile
sandArray = sandproductionrate(0.1,simLength,'cst',0.02);

% tuning (from HAC)
nmpcConfig.dumax = 0.05; % 0.01 | 0.1
nmpcConfig.umax = 2;
nmpcConfig.umin = 0.4;

% %% plotting inputs
% figure(1)
% for well = 1
%     subplot(3,1,well)
%     hold on
%     for ii = 1:nMC
%        plot(1:simLength,U_MC{well}(ii,:),'Color',[0,0,0,0.05],'Linewidth',1)
%     end
%     
%     ylim([nmpcConfig.umin nmpcConfig.umax])
%     ylabel('Q_g [L/min]','FontSize',10)
%     
%     xticks(0:50:(simLength - 1))
%     xlim([0 (simLength - 1)])
%     xlabel('time [days]','FontSize',10)
%     
%  
% end
% 
% %% plotting degradation
% figure(2)
% for m_k = 1
%     subplot(3,1,m_k)
%     hold on
%     for ii = 1:nMC
%        plot(1:simLength,erosionArray{m_k,1},'Color',[0,0,0,0.05],'Linewidth',1)
%     end
%     
%     %ylim([0 nmpcConfig.umax])
%     ylabel('Crack length [cm]','FontSize',10)
%     
%     xticks(0:50:(simLength - 1))
%     xlim([0 (simLength - 1)])
%     xlabel('time [days]','FontSize',10)
%     
% 
% end

% %% plotting inputs
% figure(1)
% for well = 1
%     subplot(3,1,well)
%     hold on
%     for ii = 1:nMC
%        plot(1:simLength,U_MC{well}(ii,:),'Color',[0,0,0,0.05],'Linewidth',1)
%     end
%     
%     ylim([nmpcConfig.umin nmpcConfig.umax])
%     ylabel('Q_g [L/min]','FontSize',10)
%     
%     xticks(0:50:(simLength - 1))
%     xlim([0 (simLength - 1)])
%     xlabel('time [days]','FontSize',10)
%     
%  
% end
% 
%% plotting ESCAPE
figure(3)
    
subplot(3,2,1)
    histogram(eps_MC(:,1),50)
    
    ylabel('Count','FontSize',10)
    xlabel('Degradation at t_0 [mm]','FontSize',10)
    title('Histogram of d_0')
    
subplot(3,2,2)
    histogram(t_MC(:,1),50)
    
    ylabel('Count','FontSize',10)
    xlabel('t_{d} [day]','FontSize',10)
    title('Histogram of t_{d}')
    
subplot(3,2,[3 4])
    hold on
    for ii = 1:nMC
       plot(1:simLength,U_MC{1}(ii,:),'Color',[0 0.4470 0.7410,0.05],'Linewidth',1,'handlevisibility','off')
    end
    
    xline(100,'-','color',[0 0.4470 0.7410],'linewidth',1);
    xline(200,'-.','color',[0 0.4470 0.7410],'linewidth',1);
    
    yline(nmpcConfig.umin,':','linewidth',2);
    yline(nmpcConfig.umax,':','linewidth',2,'handlevisibility','off')
    
    box on
    
    ylim([nmpcConfig.umin - 0.5, nmpcConfig.umax + 0.5])
    ylabel('Q_{gl} [L/min]','FontSize',10)
    
    xticks(0:50:(simLength - 1))
    xlim([0 (simLength - 1)])
    xlh = xlabel('time [days]','FontSize',10);
    
    xlh.Position = [281.7766402701628,-0.3];
    
    legend('t = 100 [days]','t = 200 [days]','Q_{gl}^{bounds}')
    
    title('Manipulated variable (gas lift injection) sequence U')

subplot(3,2,5)
    hold on
    h1 = histogram(erosionArray{1,1}(:,100),50);
    h1.FaceColor = [0.8500 0.3250 0.0980];
    h2 = histogram(erosionArray{2,1}(:,100),50);
    h2.FaceColor = [0.9290 0.6940 0.1250];
    h3 = histogram(erosionArray{3,1}(:,100),50);
    h3.FaceColor = [0.4660 0.6740 0.1880];
    
    box on
    
%     xline(mean(erosionArray{1,1}(:,200)),':','color',[0.8500 0.3250 0.0980],'linewidth',2);
%     xline(mean(erosionArray{2,1}(:,200)),':','color',[0.9290 0.6940 0.1250],'linewidth',2);
%     xline(mean(erosionArray{3,1}(:,200)),':','color',[0.4660 0.6740 0.1880],'linewidth',2);
    
    ylim([0, 40])
    
    legend('True','LR','NNR')
    
    ylabel('Count','FontSize',10)
    xlabel('Degradation at d(100) [mm]','FontSize',10)
    title('Histogram of d(t = 100 days)')
    
subplot(3,2,6)
    hold on
    h1 = histogram(erosionArray{1,1}(:,200),50);
    h1.FaceColor = [0.8500 0.3250 0.0980];
    h2 = histogram(erosionArray{2,1}(:,200),50);
    h2.FaceColor = [0.9290 0.6940 0.1250];
    h3 = histogram(erosionArray{3,1}(:,200),50);
    h3.FaceColor = [0.4660 0.6740 0.1880];
    
    box on
    
%     xline(mean(erosionArray{1,1}(:,200)),':','color',[0.8500 0.3250 0.0980],'linewidth',2);
%     xline(mean(erosionArray{2,1}(:,200)),':','color',[0.9290 0.6940 0.1250],'linewidth',2);
%     xline(mean(erosionArray{3,1}(:,200)),':','color',[0.4660 0.6740 0.1880],'linewidth',2);
    
    ylim([0, 40])
    
    legend('True','LR','NNR')
    
    ylabel('Count','FontSize',10)
    xlabel('Degradation at d(200) [mm]','FontSize',10)
    title('Histogram of d(t = 200 days)')
    

