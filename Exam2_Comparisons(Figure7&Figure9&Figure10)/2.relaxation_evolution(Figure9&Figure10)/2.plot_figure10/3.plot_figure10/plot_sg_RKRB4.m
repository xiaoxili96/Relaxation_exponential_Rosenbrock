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

load('RK4_RB4_tau=0.25_alpha=100.mat');
axes(p(1));    % 获取当前figure的信息
plot(tmesh_RB4_4(1:end-1),GAMMA_RB4_4,'-bo','LineWidth',2,'MarkerSize',8);  hold on;
plot(tmesh_RK4_4(1:end-1),GAMMA_RK4_4,'-.rx','LineWidth',2,'MarkerSize',8);
box on;  set(gca,'linewidth',2,'FontSize',24);
set(gca,'XLim',[0 2]);  set(gca,'YLim',[0 2]);
xlabel('$t$','interpreter','latex','FontSize',30);
ylabel('$\gamma_n$','interpreter','latex','FontSize',30);
legend('RER(4,3) with $\tau=1/4$','RRK(4,4) with $\tau=1/4$');
set(legend,'interpreter','latex','location','northeast','FontSize',24);