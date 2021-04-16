addpath('analysis')
% tmp = load('respMat_080219_s01.mat');
tmp = load('\\sccn.ucsd.edu\projects\Shawn\2019_VEP\respMat_1455_NSFVEP_927');
respMat = tmp.respMat{2};

vep_behavior(respMat)

% task-related component analysis
temp = load('\\sccn.ucsd.edu\projects\Shawn\2019_VEP\1455_NSFVEP_927.mat'); % epochedEEG
epochedEEG = temp.epochedEEG;
clear temp;

mean_acc = cell(5,3);
mu_ci = cell(5,3);
for STIM = 1:5
    for CONT = 1:3
        try
        [mean_acc{STIM,CONT}, mu_ci{STIM,CONT}] = vep_trca(epochedEEG,STIM,CONT);
        end
    end
end

figure,
for CONT = 1:3
    subplot(1,3,CONT); hold on,
    for STIM = 1:5
        try
        plot([0.2:0.2:3], mean_acc{STIM,CONT},'linewidth',2)
        end
    end
    ylim([0 100]); xlabel('Data Length (sec)'); ylabel('Accuracy');
    hold off
end
legend('HF-MSEQ', 'LF-MSEQ', 'LF-SSVEP','HF-SSVEP','LF-SSVEP-IMG')
set(gcf,'position',[50 50 1200 400]);


figure
hold on
ha = bar(1:nStim,mean_acc);
% For each set of bars, find the centers of the bars, and write error bars
pause(0.1); %pause allows the figure to be created
for ib = 1:numel(ha)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = ha(ib).XData+ha(ib).XOffset;
    errorbar(xData,mean_acc(:,ib),mean_acc(:,ib)-mu_ci(:,ib,1),mu_ci(:,ib,2)-mean_acc(:,ib),'k.')
end
set(gca,'XTick',[1 2 3 4 5],'XTickLabel',{'FMC', 'MSEQ', 'SSVEP', 'FMC-IMG'},'FontSize',12);
xlabel('Stimuli Types'); ylabel('Cross validation accuracy');
legend('2','8','16');

% figure, plot(mean_acc','linewidth',2);
% legend('FMC','MSEQ','SSVEP','FMC-IMG');
% xlabel('Contrast'); set(gca,'XTick',[1,2,3],'XTickLabel',{'2','4','8'},'fontsize',14);
% ylabel('Cross Validation Accuracy'); ylim([0 100]);


% reaction time - indicator of learning and fatigue?
figure, plot(respMat.rt)