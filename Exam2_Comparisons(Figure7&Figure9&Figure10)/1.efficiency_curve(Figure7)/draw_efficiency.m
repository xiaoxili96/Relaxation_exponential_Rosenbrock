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
load('data.mat');

axes(p(1));    % 获取当前figure的信息
loglog(TIME_AVF,ERR_AVF,'-y+','LineWidth',2,'MarkerSize',8);  hold on;
loglog(TIME_EAVF,ERR_EAVF,'-gx','LineWidth',2,'MarkerSize',8);
loglog(TIME_csRK4(1:5),ERR_csRK4(1:5),'-r*','LineWidth',2,'MarkerSize',8);
loglog(TIME_rexprb2s1,ERR_rexprb2s1,'-bD','LineWidth',2,'MarkerSize',8);
loglog(TIME_rexprb4s3,ERR_rexprb4s3,'-ko','LineWidth',2,'MarkerSize',8);
box on;  set(gca,'linewidth',2,'FontSize',24);
set(gca,'XLim',[0.1 1000]);  set(gca,'YLim',[10^(-9) 0.1]);
xlabel('$CPU$ $time$','interpreter','latex','FontSize',30);
ylabel('$Error$','interpreter','latex','FontSize',30);
legend('AVF','EAVF','csRK4','RER(2,1)','RER(4,3)');
set(legend,'interpreter','latex','location','southwest','FontSize',24);