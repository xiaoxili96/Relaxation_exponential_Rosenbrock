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
load('data_AVF.mat');
load('data_EAVF.mat');
load('data_csRK4.mat');
load('data_rexprb2s1.mat');
load('data_rexprb4s3.mat');

axes(p(1));    
loglog(TIME_AVF,ERR_AVF,'-y+','LineWidth',2,'MarkerSize',8);  hold on;
loglog(TIME_EAVF,ERR_EAVF,'-gx','LineWidth',2,'MarkerSize',8);
loglog(TIME_csRK4,ERR_csRK4,'-r*','LineWidth',2,'MarkerSize',8);
loglog(TIME_rexprb2s1,ERR_rexprb2s1,'-bD','LineWidth',2,'MarkerSize',8);
loglog(TIME_rexprb4s3,ERR_rexprb4s3,'-ko','LineWidth',2,'MarkerSize',8);
box on;  set(gca,'linewidth',2,'FontSize',24);
xlabel('$CPU$ $time$','interpreter','latex','FontSize',30);
ylabel('$Error$','interpreter','latex','FontSize',30);
legend('AVF','EAVF','csRK4','RER(2,1)','RER(4,3)');
set(legend,'interpreter','latex','location','northeast','FontSize',24); 