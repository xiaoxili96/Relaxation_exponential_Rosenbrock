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
load('rexprb2s1.mat');
load('rexprb3s2.mat');
load('rexprb4s3.mat');


axes(p(1));    % 获取当前figure的信息
semilogy(tmeshr2s1(1:end-1),abs(Energyr2s1(1:end-1)-Energyr2s1(1))/abs(Energyr2s1(1)),'-','LineWidth',2); hold on;
semilogy(tmeshr3s2(1:end-1),abs(Energyr3s2(1:end-1)-Energyr3s2(1))/abs(Energyr3s2(1)),'-','LineWidth',2);
semilogy(tmeshr4s3(1:end-1),abs(Energyr4s3(1:end-1)-Energyr4s3(1))/abs(Energyr4s3(1)),'-','LineWidth',2);
box on;  set(gca,'linewidth',2,'FontSize',24);
set(gca,'XLim',[0 300]);  % set(gca,'YLim',[10^(-17) 10^(-14)]);
xlabel('$t$','interpreter','latex','FontSize',30);
ylabel('$D^n_\gamma$','interpreter','latex','FontSize',30);
legend('RER(2,1) with $\tau=1/16$','RER(3,2) with $\tau=1/8$','RER(4,3) with $\tau=1/4$');  
set(legend,'interpreter','latex','location','northwest','FontSize',24); 