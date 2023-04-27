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
load('exprb2s1_100.mat');
tmesh_2s1=tmesh;  Energy_2s1=Energy;
load('exprb3s2_20.mat');
tmesh_3s2=tmesh;  Energy_3s2=Energy;
load('exprb4s3_20.mat');
tmesh_4s3=tmesh;  Energy_4s3=Energy;
load('rexprb2s1_10.mat');
tmesh_r2s1=tmesh;  Energy_r2s1=Energy;
load('rexprb3s2_2.mat');
tmesh_r3s2=tmesh;  Energy_r3s2=Energy;
load('rexprb4s3_2.mat');
tmesh_r4s3=tmesh;  Energy_r4s3=Energy;

axes(p(1));    
semilogy(tmesh_2s1,abs(Energy_2s1-Energy_2s1(1))/abs(Energy_2s1(1)),'r','LineWidth',2);  hold on;
semilogy(tmesh_r2s1,abs(Energy_r2s1-Energy_r2s1(1))/abs(Energy_r2s1(1)),'r-.','LineWidth',2);
semilogy(tmesh_3s2,abs(Energy_3s2-Energy_3s2(1))/abs(Energy_3s2(1)),'g','LineWidth',2);
semilogy(tmesh_r3s2,abs(Energy_r3s2-Energy_r3s2(1))/abs(Energy_r3s2(1)),'g-.','LineWidth',2);
semilogy(tmesh_4s3,abs(Energy_4s3-Energy_4s3(1))/abs(Energy_4s3(1)),'b','LineWidth',2);
semilogy(tmesh_r4s3,abs(Energy_r4s3-Energy_r4s3(1))/abs(Energy_r4s3(1)),'b-.','LineWidth',2);
box on;  set(gca,'linewidth',2,'FontSize',24);
set(gca,'XLim',[0 200]);  set(gca,'YLim',[10^(-16) 0.001]);
xlabel('$t$','interpreter','latex','FontSize',30);
ylabel('$D^n$ or $D^n_\gamma$','interpreter','latex','FontSize',30);
legend('ER(2,1): $\tau=1/100$, $t_c=3179s$','RER(2,1): $\tau=1/10$, $t_c=354s$','ER(3,2): $\tau=1/20$, $t_c=1308s$','RER(3,2): $\tau=1/2$, $t_c=194$','ER(4,3): $\tau=1/20$, $t_c=1970s$','RER(4,3): $\tau=1/2$, $t_c=255s$');
set(legend,'interpreter','latex','location','east','FontSize',24); 