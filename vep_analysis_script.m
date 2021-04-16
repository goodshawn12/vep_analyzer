%addpath('analysis')
tmp = load('respMat_1455_NSFVEP_s3.mat');
respMat = tmp.respMat{2};

figure();
vep_behavior_img(respMat)
saveas(gcf,[pwd '/VEP experiment 2.0/analysis/S1/Response.jpg']);

% task-related component analysis
temp = load('NSFVEP_s3.mat'); % epochedEEG
epochedEEG = temp.epochedEEG;


clear temp;

epochedEEG_truc = epochedEEG(1:4,:,:,:);
mean_acc = zeros(4,3);
mu_ci = zeros(4,3,2);
for STIM = 1:4
    for CONT = 1:3
        [mean_acc(STIM,CONT), mu_ci(STIM,CONT,:)] = vep_trca(epochedEEG_truc,STIM,CONT);
    end
end


figure
hold on
ha = bar(1:STIM,mean_acc);
% For each set of bars, find the centers of the bars, and write error bars
pause(0.1); %pause allows the figure to be created
for ib = 1:numel(ha)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = ha(ib).XData+ha(ib).XOffset;
    errorbar(xData,mean_acc(:,ib),mean_acc(:,ib)-mu_ci(:,ib,1),mu_ci(:,ib,2)-mean_acc(:,ib),'k.')
end
set(gca,'XTick',[1 2 3 4 5],'XTickLabel',{'FMC', 'MSEQ', 'LF-SSVEP', 'HF-SSVEP'},'FontSize',12);
xlabel('Stimuli Types'); ylabel('Cross validation accuracy');
legend('2','8','16');
saveas(gcf,[pwd '/VEP experiment 2.0/analysis/S3/Accuracy.jpg']);

%%
figure, plot(mean_acc','linewidth',2);
 legend('FMC','MSEQ','SSVEP','FMC-IMG');
 xlabel('Contrast'); set(gca,'XTick',[1,2,3],'XTickLabel',{'2','4','8'},'fontsize',14);
 ylabel('Cross Validation Accuracy'); ylim([0 100]);


% reaction time - indicator of learning and fatigue?
%figure, plot(respMat.rt)



%%

mean_acc_var = cell(4,3);
mu_ci_var = cell(4,3);
cont = ['2','8','16'];
for CONT = 1:3
    cont_i = cont(CONT);
    fig = figure();
    hold on
    for STIM = 1:4
        [mean_acc_var{STIM,CONT}, mu_ci_var{STIM,CONT}] = vep_trca_var(epochedEEG_truc,STIM,CONT);
    end
    hold off
    saveas(gcf,[pwd '/VEP experiment 2.0/analysis/S3/' cont_i '.jpg']);
end

%%

data_length = [0.2:0.1:3.4];
mean_acc_var = cell(4,3,10);
mu_ci_var = cell(4,3,10);
for SUB = 1:10
    name = sprintf('NSFVEP_s%d.mat', SUB);
    temp = load(name); % epochedEEG
    epochedEEG = temp.epochedEEG;
    clear temp;
    epochedEEG_truc = epochedEEG(1:4,:,:,:); 
    
    for CONT = 1:3
        for STIM = 1:4
        [mean_acc_var{STIM,CONT,SUB}, mu_ci_var{STIM,CONT,SUB}] = vep_trca_var(epochedEEG_truc,STIM,CONT);
        end
    end
end

crosssub_mean_var = cell(3);
crosssub_sd_var = cell(3);

for CONT = 1:3
    crosssub_mean_stim = zeros(4,length(data_length));
    crosssub_sd_stim = zeros(4,length(data_length));
    for STIM = 1:4
        temp = zeros(length(data_length),10);
        for SUB = 1:10
            temp(:,SUB) = mean_acc_var{STIM,CONT,SUB};
        end
        crosssub_mean_stim(STIM,:) = mean(temp,2);
        crosssub_sd_stim(STIM,:) = std(temp,[],2);
    end
    crosssub_mean_var{CONT} = crosssub_mean_stim;
    crosssub_sd_var{CONT} = crosssub_sd_stim;
end

cont = [2,8,16];
for CONT_FIG = 1:3
    fig = figure();
    crosssub_mean_stim_var = crosssub_mean_var{CONT_FIG};
    crosssub_sd_stim_var = crosssub_sd_var{CONT_FIG};
    hold on
    for STIM_FIG = 1:4
       errorbar(data_length, crosssub_mean_stim_var(STIM_FIG, :), crosssub_sd_stim_var(STIM_FIG, :), 'linewidth',2);
    end
    hold off
     xlabel('Training data legnth (sec)'); ylabel('Cross validation accuracy (''%)');
     legend('FMC','MSEQ','LF_SSVEP','HF_SSVEP');
     legend('Location', 'northwest');
     set(gca,'fontsize',12)
     num = sprintf('/VEP experiment 2.0/analysis/%d.jpg', cont(CONT_FIG));
     saveas(gcf,[pwd num]);
end

