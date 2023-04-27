%% INPUT TIGHTPLOT PARAMETERS
TightPlot.ColumeNumber = 1;     % ��ͼ����
TightPlot.RowNumber = 1;    % ��ͼ����
TightPlot.GapW = 0.05;  % ��ͼ֮������Ҽ��
TightPlot.GapH = 0.05;   % ��ͼ֮������¼��
TightPlot.MarginsLower = 0.13;   % ��ͼ��ͼƬ�·��ļ��
TightPlot.MarginsUpper = 0.11;  % ��ͼ��ͼƬ�Ϸ��ļ��
TightPlot.MarginsLeft = 0.24;   % ��ͼ��ͼƬ�󷽵ļ��
TightPlot.MarginsRight = 0.38;  % ��ͼ��ͼƬ�ҷ��ļ��

%% PLOT
figure(1);  % ����Figure
p = tight_subplot(TightPlot.ColumeNumber,TightPlot.RowNumber,...
    [TightPlot.GapH TightPlot.GapW],...
    [TightPlot.MarginsLower TightPlot.MarginsUpper],...
    [TightPlot.MarginsLeft TightPlot.MarginsRight]);    % �������ò�����һ���Ѿ�������

load('RK4_RB4_tau=0.25_alpha=100.mat');
axes(p(1));    % ��ȡ��ǰfigure����Ϣ
plot(tmesh_RB4_4(1:end-1),GAMMA_RB4_4,'-bo','LineWidth',2,'MarkerSize',8);  hold on;
plot(tmesh_RK4_4(1:end-1),GAMMA_RK4_4,'-.rx','LineWidth',2,'MarkerSize',8);
box on;  set(gca,'linewidth',2,'FontSize',24);
set(gca,'XLim',[0 2]);  set(gca,'YLim',[0 2]);
xlabel('$t$','interpreter','latex','FontSize',30);
ylabel('$\gamma_n$','interpreter','latex','FontSize',30);
legend('RER(4,3) with $\tau=1/4$','RRK(4,4) with $\tau=1/4$');
set(legend,'interpreter','latex','location','northeast','FontSize',24);