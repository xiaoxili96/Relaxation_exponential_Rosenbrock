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
load('data.mat');

N=40;  h=2*pi/N;  xmesh=0:h:2*pi-h;

axes(p(1));    
[XX,YY]=meshgrid(tmesh,[xmesh 2*pi]);  surf(XX,YY,[U_save;U_save(1,:)]);  shading interp;  colormap Jet;
box on;  set(gca,'linewidth',2,'FontSize',24);
set(gca,'XLim',[0 2])
set(gca,'YLim',[0 2*pi])
set(gca,'ZLim',[-4 4])
xlabel('$t$','interpreter','latex','FontSize',30);
ylabel('$x$','interpreter','latex','FontSize',30);
zlabel('$U$','interpreter','latex','FontSize',30);