%% INPUT TIGHTPLOT PARAMETERS
TightPlot.ColumeNumber = 1;     % 子图行数
TightPlot.RowNumber = 1;    % 子图列数
TightPlot.GapW = 0.05;  % 子图之间的左右间距
TightPlot.GapH = 0.05;   % 子图之间的上下间距
TightPlot.MarginsLower = 0.13;   % 子图与图片下方的间距
TightPlot.MarginsUpper = 0.11;  % 子图与图片上方的间距
TightPlot.MarginsLeft = 0.24;   % 子图与图片左方的间距
TightPlot.MarginsRight = 0.38;  % 子图与图片右方的间距

%% PLOT
figure(1);  % 声明Figure
p = tight_subplot(TightPlot.ColumeNumber,TightPlot.RowNumber,...
    [TightPlot.GapH TightPlot.GapW],...
    [TightPlot.MarginsLower TightPlot.MarginsUpper],...
    [TightPlot.MarginsLeft TightPlot.MarginsRight]);    % 具体设置参数上一节已经输入了

N=40;  h=2*pi/N;  xmesh=0:h:2*pi-h;

axes(p(1));    % 获取当前figure的信息
[XX,YY]=meshgrid(tmesh,[xmesh 2*pi]);  surf(XX,YY,[U_save;U_save(1,:)]);  shading interp;  colormap Jet;
box on;  set(gca,'linewidth',2,'FontSize',24);
set(gca,'XLim',[8 10])
set(gca,'YLim',[0 2*pi])
set(gca,'ZLim',[-4 4])
xlabel('$t$','interpreter','latex','FontSize',30);
ylabel('$x$','interpreter','latex','FontSize',30);
zlabel('$U$','interpreter','latex','FontSize',30);