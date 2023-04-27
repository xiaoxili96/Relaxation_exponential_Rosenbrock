%% INPUT TIGHTPLOT PARAMETERS
TightPlot.ColumeNumber = 2;     
TightPlot.RowNumber = 3;    
TightPlot.GapW = 0.1;  
TightPlot.GapH = 0.1;   
TightPlot.MarginsLower = 0.06;  
TightPlot.MarginsUpper = 0.04; 
TightPlot.MarginsLeft = 0.08;   
TightPlot.MarginsRight = 0.158;  

%% PLOT
figure(1);  
p = tight_subplot(TightPlot.ColumeNumber,TightPlot.RowNumber,...
    [TightPlot.GapH TightPlot.GapW],...
    [TightPlot.MarginsLower TightPlot.MarginsUpper],...
    [TightPlot.MarginsLeft TightPlot.MarginsRight]);    
load('rexprb4s3.mat');

N=40;  left=-7;  right=7;  bottom=-7; top=7;  h=(right-left)/N;  
xmesh=left+0.5*h:h:right-0.5*h;  ymesh=xmesh;  [XMESH,YMESH]=meshgrid(xmesh,ymesh); 

axes(p(1));    
Un=U_save(:,1);
surf(XMESH,YMESH,reshape(Un,N,N));  shading interp;  colormap Jet;  axis([-7,7,-7,7,-0.5,1]);  drawnow;
box on;  set(gca,'linewidth',2,'FontSize',16);
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
zlabel('$\sin(U/2)$','interpreter','latex');
title('$t=0$','interpreter','latex');
axes(p(2));    
Un=U_save(:,27);
surf(XMESH,YMESH,reshape(Un,N,N));  shading interp;  colormap Jet;  axis([-7,7,-7,7,-0.5,1]);  drawnow;
box on;  set(gca,'linewidth',2,'FontSize',16);
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
zlabel('$\sin(U/2)$','interpreter','latex');
title('$t=2.6$','interpreter','latex');
axes(p(3));    
Un=U_save(:,57);
surf(XMESH,YMESH,reshape(Un,N,N));  shading interp;  colormap Jet;  axis([-7,7,-7,7,-1.5,1]);  drawnow;
box on;  set(gca,'linewidth',2,'FontSize',16);
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
zlabel('$\sin(U/2)$','interpreter','latex');
title('$t=5.6$','interpreter','latex');
axes(p(4));    
Un=U_save(:,85);
surf(XMESH,YMESH,reshape(Un,N,N));  shading interp;  colormap Jet;  axis([-7,7,-7,7,-1,0.5]);  drawnow;
box on;  set(gca,'linewidth',2,'FontSize',16);
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
zlabel('$\sin(U/2)$','interpreter','latex');
title('$t=8.4$','interpreter','latex');
axes(p(5));    
Un=U_save(:,113);
surf(XMESH,YMESH,reshape(Un,N,N));  shading interp;  colormap Jet;  axis([-7,7,-7,7,-1.5,1]);  drawnow;
box on;  set(gca,'linewidth',2,'FontSize',16);
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
zlabel('$\sin(U/2)$','interpreter','latex');
title('$t=11.2$','interpreter','latex');
axes(p(6));    
Un=U_save(:,127);
surf(XMESH,YMESH,reshape(Un,N,N));  shading interp;  colormap Jet;  axis([-7,7,-7,7,-1.5,1]);  drawnow;
box on;  set(gca,'linewidth',2,'FontSize',16);
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
zlabel('$\sin(U/2)$','interpreter','latex');
title('$t=12.6$','interpreter','latex');