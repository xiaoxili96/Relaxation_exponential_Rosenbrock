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
load('rexprb2s1_50.mat');
tmesh_r2s1=tmesh;  Energy_r2s1=Energy;
load('rexprb3s2_45.mat');
tmesh_r3s2=tmesh;  Energy_r3s2=Energy;
load('rexprb4s3_20.mat');
tmesh_r4s3=tmesh;  Energy_r4s3=Energy;

axes(p(1));    
semilogy(tmesh_r2s1(1:end-1),abs(Energy_r2s1(1:end-1)-Energy_r2s1(1))/abs(Energy_r2s1(1)),'-','LineWidth',2); hold on;
semilogy(tmesh_r3s2(1:end-1),abs(Energy_r3s2(1:end-1)-Energy_r3s2(1))/abs(Energy_r3s2(1)),'-','LineWidth',2);
semilogy(tmesh_r4s3(1:end-1),abs(Energy_r4s3(1:end-1)-Energy_r4s3(1))/abs(Energy_r4s3(1)),'-','LineWidth',2);
box on;  set(gca,'linewidth',2,'FontSize',24);
set(gca,'XLim',[0 200]);  % set(gca,'YLim',[10^(-16) 10^(-13)]);
xlabel('$t$','interpreter','latex','FontSize',30);
ylabel('$D^n_\gamma$','interpreter','latex','FontSize',30);
legend('RER(2,1) with $\tau=1/50$','RER(3,2) with $\tau=1/45$','RER(4,3) with $\tau=1/20$');
set(legend,'interpreter','latex','location','northwest','FontSize',24);