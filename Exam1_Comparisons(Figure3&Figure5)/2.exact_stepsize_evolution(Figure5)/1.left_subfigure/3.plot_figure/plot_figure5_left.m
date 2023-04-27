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
load('data_rexprb2s1.mat');
load('data_rrk2s2.mat');

axes(p(1));    % 获取当前figure的信息 
semilogy(tmesh_rexprb2s1_2(1:end-1),exact_stepsize_rexprb2s1_2,'-rs','LineWidth',2,'MarkerSize',8);  hold on; 
semilogy(tmesh_rexprb2s1_4(1:end-1),exact_stepsize_rexprb2s1_4,'-kD','LineWidth',2,'MarkerSize',8);
semilogy(tmesh_rexprb2s1_8(1:end-1),exact_stepsize_rexprb2s1_8,'-bo','LineWidth',2,'MarkerSize',8);
semilogy(tmesh_rrk2s2_2(1:end-1),exact_stepsize_rrk2s2_2,'-.y*','LineWidth',2,'MarkerSize',8);  
semilogy(tmesh_rrk2s2_4(1:end-1),exact_stepsize_rrk2s2_4,'-.m+','LineWidth',2,'MarkerSize',8);
semilogy(tmesh_rrk2s2_8(1:end-1),exact_stepsize_rrk2s2_8,'-.gx','LineWidth',2,'MarkerSize',8);
box on;  set(gca,'linewidth',2,'FontSize',24);
set(gca,'XLim',[0 2]);  set(gca,'YLim',[0.00001 1]);
xlabel('$t$','interpreter','latex','FontSize',30);
ylabel('$\gamma_n \tau$','interpreter','latex','FontSize',30);
legend('RER(2,1) with $\tau=1/2$','RER(2,1) with $\tau=1/4$','RER(2,1) with $\tau=1/8$','RRK(2,2) with $\tau=1/2$','RRK(2,2) with $\tau=1/4$','RRK(2,2) with $\tau=1/8$')
set(legend,'interpreter','latex','location','southeast','FontSize',24);