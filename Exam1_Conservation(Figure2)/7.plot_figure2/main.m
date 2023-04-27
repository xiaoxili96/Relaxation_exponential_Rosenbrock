%% INPUT TIGHTPLOT PARAMETERS
TightPlot.ColumeNumber = 1;    
TightPlot.RowNumber = 1;    
TightPlot.GapW = 0.05;  
TightPlot.GapH = 0.05;   
TightPlot.MarginsLower = 0.13;  
TightPlot.MarginsUpper = 0.11;  
TightPlot.MarginsLeft = 0.24;  
TightPlot.MarginsRight = 0.38; 

%% PLOT
figure(1);  
p = tight_subplot(TightPlot.ColumeNumber,TightPlot.RowNumber,...
    [TightPlot.GapH TightPlot.GapW],...
    [TightPlot.MarginsLower TightPlot.MarginsUpper],...
    [TightPlot.MarginsLeft TightPlot.MarginsRight]);   
load('exprb2s1.mat');
load('exprb3s2.mat');
load('exprb4s3.mat');
load('rexprb2s1.mat');
load('rexprb3s2.mat');
load('rexprb4s3.mat');


axes(p(1));    % 获取当前figure的信息
semilogy(tmesh2s1,abs(Energy2s1-Energy2s1(1))/abs(Energy2s1(1)),'r','LineWidth',2);  hold on;
semilogy(tmeshr2s1(1:end-1),abs(Energyr2s1(1:end-1)-Energyr2s1(1))/abs(Energyr2s1(1)),'r-.','LineWidth',2);
semilogy(tmesh3s2,abs(Energy3s2-Energy3s2(1))/abs(Energy3s2(1)),'g','LineWidth',2);
semilogy(tmeshr3s2(1:end-1),abs(Energyr3s2(1:end-1)-Energyr3s2(1))/abs(Energyr3s2(1)),'g-.','LineWidth',2);
semilogy(tmesh4s3,abs(Energy4s3-Energy4s3(1))/abs(Energy4s3(1)),'b','LineWidth',2);
semilogy(tmeshr4s3(1:end-1),abs(Energyr4s3(1:end-1)-Energyr4s3(1))/abs(Energyr4s3(1)),'b-.','LineWidth',2);
box on;  set(gca,'linewidth',2,'FontSize',24);
set(gca,'XLim',[0 200]);  set(gca,'YLim',[10^(-16) 1]);
xlabel('$t$','interpreter','latex','FontSize',30);
ylabel('$D^n$ or $D^n_\gamma$','interpreter','latex','FontSize',30);
legend('ER(2,1): $\tau=1/500$, $t_c=304s$','RER(2,1): $\tau=1/50$, $t_c=44s$','ER(3,2): $\tau=1/300$, $t_c=324s$','RER(3,2): $\tau=1/30$, $t_c=56s$','ER(4,3): $\tau=1/200$, $t_c=322s$','RER(4,3): $\tau=1/20$, $t_c=50s$');  
set(legend,'interpreter','latex','location','east','FontSize',24); 