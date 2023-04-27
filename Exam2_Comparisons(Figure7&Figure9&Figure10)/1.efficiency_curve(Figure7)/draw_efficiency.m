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
load('data.mat');

axes(p(1));    % ��ȡ��ǰfigure����Ϣ
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