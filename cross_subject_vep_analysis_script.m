% ------------------------------------------------------------------------
%               Cross Subject Accuracy Analysis Script
% ------------------------------------------------------------------------
data_length = [0.1:0.1:3.0];
CH_label = [2,3,4,9,16,25,32];
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
    crosssub_mean_stim = zeros(4,length(CH_label));
    crosssub_sd_stim = zeros(4,length(CH_label));
    for STIM = 1:4  %extract and determine the mean and standard deviation for each data point
        temp = zeros(length(CH_label),10);
        for SUB = 1:10
            temp(:,SUB) = mean_acc_var{STIM,CONT,SUB};
        end
        crosssub_mean_stim(STIM,:) = mean(temp,2);
        crosssub_sd_stim(STIM,:) = std(temp,[],2)./sqrt(10);
    end
    crosssub_mean_var{CONT} = crosssub_mean_stim;   %store the data in respect to contrast
    crosssub_sd_var{CONT} = crosssub_sd_stim;
end

cont = [2,8,16];
for CONT_FIG = 1:3
    fig = figure();
    crosssub_mean_stim_var = crosssub_mean_var{CONT_FIG};
    crosssub_sd_stim_var = crosssub_sd_var{CONT_FIG};
    hold on
    for STIM_FIG = 1:4
        errorbar(CH_label, crosssub_mean_stim_var(STIM_FIG, :), crosssub_sd_stim_var(STIM_FIG, :), 'linewidth',2);
    end
    hold off
    xlabel('Number of Channels'); ylabel('Cross validation accuracy (%)');
    legend('FMC','MSEQ','LF_SSVEP','HF_SSVEP');
    legend('Location', 'southeast');
    set(gca,'fontsize',12)
    ylim([0 100]);
    num = sprintf('/VEP experiment 2.0/analysis/%d.jpg', cont(CONT_FIG));
    saveas(gcf,[pwd num]);
end

%%
% ------------------------------------------------------------------------
%               ERP Image plotting
% ------------------------------------------------------------------------
chan = [1:32];
%cont_id = 1:3;
STIM_label = {['FMC'],['MSEQ'],['LF_SSVEP'],['HF_SSVEP']};
CONT_label = {['CONT=2'], ['CONT=8'], ['CONT=16']};
erp_temp = cell(4,3,10,4);
erp_img = cell(4,3);

for sub_id = 1:10
    name = sprintf('NSFVEP_s%d.mat', sub_id);
    temp = load(name); % epochedEEG
    epochedEEG = temp.epochedEEG;
    clear temp;
    epochedEEG_truc = epochedEEG(1:4,:,:,:);
    
    for stim_id = 1:4
        for cont_id = 1:3
            for loc_id = 1:4
                for tr_id = 1:10
                    erp_temp{stim_id, cont_id,sub_id, loc_id} = [erp_temp{stim_id, cont_id, sub_id, loc_id}; mean(epochedEEG_truc{stim_id,cont_id,loc_id,tr_id}(chan,:),1)];
                end
            end
        end
    end
end

for STIM = 1:4
    for CONT = 1:3
        for LOC = 1:4
            for SUB = 1:10
                erp_img{STIM, CONT} = [erp_img{STIM, CONT}; erp_temp{STIM, CONT, SUB, LOC}];
            end
        end
    end
end

for stim_id = 1:4
    for cont_id = 1:3
        figure, imagesc(erp_img{stim_id, cont_id});
        xlim([102 1638]);
        xticks([102, 614, 1126, 1638]);
        xticklabels({'0', '1', '2', '3'});
        yticks(linspace(0, 400, 9));
        yticklabels({' ','    loc1     ', '------    ','    loc2     ','------    ','    loc3     ','------    ', '    loc4     ',' '})
        xlabel('Time (sec)'); ylabel('Location of Target');
        colorbar
        caxis([-20 20]);
        stim_name = ['/VEP experiment 2.0/Analysis/Cross_subject/' STIM_label{stim_id} '_' CONT_label{cont_id} '.jpg'];
        saveas(gcf,[pwd stim_name]);
    end
end

%% ***Previous version, not for actual use now, only saved for a record. ***
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

stat_Mat = zeros(30,4);

for STIM_id = 1:4
    for CONT_id = 1:3
        for SUB_id = 1:10
            row_index = (CONT_id-1)*10+SUB_id;
            col_index = STIM_id;
            stat_Mat(row_index,col_index) = mean_acc_var{STIM_id,CONT_id,SUB_id};
        end
    end
end

[~,~,stats] = anova2(stat_Mat);

ttest_results = cell(4,1);
cond1 = zeros(10,1);
cond2 = zeros(10,1);
for stim = 1:4
    ttest_stim_results = zeros(3,3);
    for cont = 1:3
        count = 1;
        for sub = 1:10
            cond1(sub) = mean_acc_var{stim, cont, sub};
        end
        for stim_comp = 1:4
            if stim_comp ~= stim
                for sub2 = 1:10
                    cond2(sub2) = mean_acc_var{stim_comp, cont, sub2};
                end
                [~,ttest_stim_results(count, cont)] = ttest(cond1,cond2);
                count = count + 1;
            end
        end
    end
    ttest_results{stim} = ttest_stim_results;
end
