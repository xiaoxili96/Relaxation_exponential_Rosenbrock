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
load('data_rexprb3s2.mat');
load('data_rrk3s3.mat');

axes(p(1));    
semilogy(tmesh_rexprb3s2_2(1:end-1),exact_stepsize_rexprb3s2_2,'-rs','LineWidth',2,'MarkerSize',8);  hold on;
semilogy(tmesh_rexprb3s2_4(1:end-1),exact_stepsize_rexprb3s2_4,'-kD','LineWidth',2,'MarkerSize',8);
semilogy(tmesh_rexprb3s2_8(1:end-1),exact_stepsize_rexprb3s2_8,'-bo','LineWidth',2,'MarkerSize',8);
semilogy(tmesh_rrk3s3_2(1:end-1),exact_stepsize_rrk3s3_2,'-.y*','LineWidth',2,'MarkerSize',8);   
semilogy(tmesh_rrk3s3_4(1:end-1),exact_stepsize_rrk3s3_4,'-.m+','LineWidth',2,'MarkerSize',8);
semilogy(tmesh_rrk3s3_8(1:end-1),exact_stepsize_rrk3s3_8,'-.gx','LineWidth',2,'MarkerSize',8);
box on;  set(gca,'linewidth',2,'FontSize',24);
set(gca,'XLim',[0 2]);  set(gca,'YLim',[0 1]);
xlabel('$t$','interpreter','latex','FontSize',30);
ylabel('$\gamma_n \tau$','interpreter','latex','FontSize',30);
legend('RER(3,2) with $\tau=1/2$','RER(3,2) with $\tau=1/4$','RER(3,2) with $\tau=1/8$','RRK(3,3) with $\tau=1/2$','RRK(3,3) with $\tau=1/4$','RRK(3,3) with $\tau=1/8$')
set(legend,'interpreter','latex','location','east','FontSize',24);