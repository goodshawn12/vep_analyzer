%% Main script for data analysis and figure generation for manuscript
% cd 'C:\Users\isnl\Documents\VEP\analysis'
% addpath('C:\Users\isnl\Documents\VEP\data')
cd 'D:\expVEP_PTB\analysis'
addpath('D:\expVEP_PTB\data');
addpath(genpath('lib')); 
addpath('metadata');

%% Fig 1 - Example of code sequence and EEG responses

% Example of 1-sec code sequences of four types of stimuli (pick one loc)
load('respMat_s1'); % same code for all subjects
codename = {'fmc','mseq','lfssvep','hfssvep'};
for code_i = 1:length(codename)
    eval(sprintf('code = respMat{1}.code_%s{4}(1:60);',codename{code_i}));
    code = [0.5 code];
    code_interp = interp1([1:1.0:length(code)],1.0*code,[1:0.05:length(code)],'previous');
    figure, plot(0:0.05:length(code)-1,code_interp,'linewidth',1.5); xlabel('Frames (1/60 sec)'); ylabel('Code');
    set(gca,'YLim',[-0.2 1.2],'XLim',[0 60],'YTick',[0 1], 'YTickLabel',{'-','+'},'fontsize',14);
    set(gcf,'position',[100 100 1300 250]);
    eval(sprintf('export_fig ''Fig1_%s'' -png -transparent',codename{code_i}));
end

% Frequency response of the 1-sec code sequence (optional)

%% Fig 1 - Averaged EEG across trials (or subjects) for each type at contrast 16
CONT = 3; LOC = 4; 
WINDOW = 1; 
SRATE = 500;
OFFSET = 0.3; % 200 msec baseline + estimated 100 msec delay
CH = 19; % Oz channel (largest SNR)
avg_data = zeros(10,SRATE*WINDOW,4);
for sub_id = 1:10
    fprintf('Processing subject %d...\n',sub_id);
    temp = load(sprintf('NSFVEP_s%d.mat', sub_id)); % epochedEEG
    epochedEEG = temp.epochedEEG;
    for stim_id = 1:4      
        epochedEEG_sel = squeeze(epochedEEG(stim_id,CONT,LOC,:));        
        data_all = zeros(10,SRATE*WINDOW);
        for tr_id = 1:10
            data = epochedEEG_sel{tr_id};
            time_range = (OFFSET*SRATE+1):((OFFSET+WINDOW)*SRATE);
            if length(CH) == 1
                data_all(tr_id,:) = data(CH,time_range);
            else
                data_all(tr_id,:) = mean(data(CH,time_range));
            end
        end
        avg_data(sub_id,:,stim_id) = mean(data_all);
    end
end

% plot averaged EEG across subjects for each stimuli
for stim_id = 1:4
    avg_data_all = mean(avg_data(:,:,stim_id));
    figure, plot([1:WINDOW*SRATE]/SRATE*1000, avg_data_all,'linewidth',1.5);
    xlabel('Time (msec)'); ylabel('EEG at Oz (\muV)');
    set(gca,'YLim',[-6 6],'fontsize',14);
    set(gcf,'position',[100 100 1300 400]);
    eval(sprintf('export_fig ''Fig1_eeg_%s'' -png -transparent',codename{stim_id}));
end


%% Fig 2 - subjective ratings of perception

% load ratings data (verify the code for generating this data)
load('average_resp.mat');

% switch order and compute mean and std across subjects
average_resp = average_resp([4,1:3],:,:);
avg_rating = mean(average_resp,3)';
std_rating = std(average_resp,[],3)'/sqrt(size(average_resp,3));

% bar plots
figure, hold on
ha = bar(1:3,avg_rating);
pause(0.1)
for ib = 1:numel(ha)
    xData = ha(ib).XData+ha(ib).XOffset;
    errorbar(xData,avg_rating(:,ib),std_rating(:,ib),'k.')
end
set(gca,'XTick',[1 2 3],'XTickLabel',{'2','8','16'},'FontSize',12);
xlabel('Contrast levels'); ylabel('Perceived Intensity Ratings');
legend('HF-SSVEP', 'FMSEQ', 'MSEQ', 'LF-SSVEP','location','northwest');
% export_fig 'Fig2_rating' -png -eps -transparent;

%% Fig2 - statistics
% prepare data for anova test
sub_num = 10;
row_len = 4; % Different colomns are different stimulis, which makes length of each row 4
col_len = 3; % Different colomns are different contrasts, which makes length of each column 3xsubject_number
mat_resp = zeros(row_len*sub_num,col_len);

for col = 1:col_len
    for row = 1:row_len
        for sub_n = 1:sub_num
            row_num = sub_num*(row - 1)+sub_n;
            mat_resp(row_num,col) = average_resp(row,col,sub_n);
        end
    end
end

% 3-way ANOVA (subject ID as the third variable) 
resp_col = average_resp(:);
stim_col = repmat([1:4]',3*10,1);
cont_col = repmat([ones(1,4), 2*ones(1,4), 3*ones(1,4)]',10,1);
subj_col = reshape(repmat([1:10],12,1),[],1);

[p,tbl,stats] = anovan(resp_col,{stim_col cont_col subj_col},'model','interaction',...
    'varnames',{'stim','cont','subj'});

% repeated measure 2-way ANOVA testing
% 1. ranova() in MATLAB
% 2. statcond() in EEGLAB
resp_cell = cell(4,3);
for STIM = 1:4
    for CONT = 1:3
        resp_cell{STIM,CONT} = squeeze(average_resp(STIM,CONT,:))';
    end
end

% perform a paired 2-way ANOVA
[F df pvals] = statcond(resp_cell,'paired','on','method','param','verbose','on');  
disp(pvals);
% returned results are very similar to 3-way ANOVA results

% post-hoc paired t-test for stimuli with manual multiple comparison correction
pval_stim = zeros(4,4,3);
for CONT = 1:3
    for STIM_1 = 1:4
        for STIM_2 = 1:4
            [~,pval_stim(STIM_1,STIM_2,CONT)] = ...
                ttest(average_resp(STIM_1,CONT,:), average_resp(STIM_2,CONT,:));
        end
    end
end

% % Archived
% % 2-way anova testing (not repeated measuers - inappropriate test)
% [p,tbl,stats] = anova2(mat_resp,10);
% 
% % multiple comparison 
% c_cont = multcompare(stats) % test for contrast levels
% c_stim = multcompare(stats,'Estimate','row') % test for stimuli

%% Fig 3. Decoding Accuracy using TRCA
% decoding accruacy over data length
DATA_LEN_LIST = [0.1:0.1:3.0];
CH_LIST = [12,17:21,24:26]; % 9 channels centered at Oz

% main loop - TRCA
mean_acc_var = cell(4,3,10);
mu_ci_var = cell(4,3,10);
for SUB = 1:10
    % load EEG (epochedEEG)
    fprintf('Processing Subject %d...\n',SUB);
    load(sprintf('NSFVEP_s%d.mat', SUB));
    epochedEEG_truc = epochedEEG(1:4,:,:,:);
    
    % averaged accuracy of leave-1-trial-out cross valiadtion using TRCA
    for CONT = 1:3
        for STIM = 1:4
            [mean_acc_var{STIM,CONT,SUB}, mu_ci_var{STIM,CONT,SUB}] = vep_trca_var(epochedEEG_truc,STIM,CONT,DATA_LEN_LIST,CH_LIST);
        end
    end
end

% save results
save('mean_acc_var_9ch.mat','mean_acc_var','mu_ci_var');

% average results across subjects
data_offset = 0.1; % for computing itr (in sec)
avg_mean_var = cell(1,3);
avg_sd_var = cell(1,3);
avg_itr = cell(1,3);
for CONT = 1:3
    avg_mean_stim = zeros(4,length(DATA_LEN_LIST));
    avg_sd_stim = zeros(4,length(DATA_LEN_LIST));
    avg_itr_stim = zeros(4,length(DATA_LEN_LIST));
    for STIM = 1:4  % extract and determine the mean and standard deviation for each data point
        temp = zeros(length(DATA_LEN_LIST),10);
        temp_itr = zeros(length(DATA_LEN_LIST),10);
        for SUB = 1:10
            temp(:,SUB) = mean_acc_var{STIM,CONT,SUB};
            for LEN = 1:length(DATA_LEN_LIST)
                temp_itr(LEN,SUB) = itr(4, mean_acc_var{STIM,CONT,SUB}(LEN)/100, DATA_LEN_LIST(LEN)+data_offset) ;
            end
        end
        avg_itr_stim(STIM,:) = mean(temp_itr,2);
        avg_mean_stim(STIM,:) = mean(temp,2);
        avg_sd_stim(STIM,:) = std(temp,[],2)./sqrt(10);
    end
    avg_itr{CONT} = avg_itr_stim([4,1:3],:);
    avg_mean_var{CONT} = avg_mean_stim([4,1:3],:);   % switch order of stimuli to match with Fig. 2
    avg_sd_var{CONT} = avg_sd_stim([4,1:3],:);
end

cont = [2,8,16];
for CONT_FIG = 1:3
    fig = figure(); hold on,
    for STIM_FIG = 1:4
        plot(DATA_LEN_LIST, avg_mean_var{CONT_FIG}(STIM_FIG,:),'linewidth',2);
    end   
    xlabel('Data length (sec)'); ylabel('Accuracy (%)');
    legend('HF-SSVEP', 'FMSEQ', 'MSEQ', 'LF-SSVEP','location','northeastoutside');
    set(gca,'fontsize',14)
    set(gcf,'position',[50, 50, 750, 500])
    ylim([0 100]);
%     export_fig(sprintf('Fig3_rating_cont%d',cont(CONT_FIG)),'-png','-transparent');close
end


%% Fig 3. Decoding Accuracy - Statistical Analysis


%% Fig 4. Information Transfer Rate (ITR)
cont = [2,8,16];
for CONT_FIG = 1:3
    fig = figure(); hold on,
    for STIM_FIG = 1:4
        plot(DATA_LEN_LIST, avg_itr{CONT_FIG}(STIM_FIG,:),'linewidth',2);
    end   
    xlabel('Data length (sec)'); ylabel('Information Transfer Rate (bits/sec)');
    legend('HF-SSVEP', 'FMSEQ', 'MSEQ', 'LF-SSVEP','location','northeastoutside');
    set(gca,'fontsize',14)
    set(gcf,'position',[50, 50, 750, 500])
    ylim([0 400])
    export_fig(sprintf('Fig4_ITR_cont%d_fixed',cont(CONT_FIG)),'-png','-transparent');close
end


%% Fig 5. Scatter plot of decoding accuracy & perception
acc = load('mean_acc_var_9ch.mat');
resp = load('average_resp.mat');

% z-scored accuracy within each subject
select_data_len = 1.0;
data_len = [0.1:0.1:3];
index_len = find(data_len == select_data_len);

mean_acc_sel_len = zeros(4,3,10);
for stim = 1:4
    for cont = 1:3
        for subj = 1:10
            mean_acc_sel_len(stim,cont,subj) = acc.mean_acc_var{stim,cont,subj}(index_len);
        end
    end
end

%% Fig 5-1: scatter plot to raw or z-scored accuracy / ratings
IS_ZSCORE = 1;
IS_LINEFIT = 1;
cont_list = [1,1,1,1,2,2,2,2,3,3,3,3];
stim_list = [1,2,3,4,1,2,3,4,1,2,3,4];
cont_color = {[0         0.4470    0.7410], ...
              [0.8500    0.3250    0.0980], ...
              [0.9290    0.6940    0.1250], ...
              [0.4940    0.1840    0.5560]};
stim_symbol = {'^','o','+','x'};

% scatter plot
acc_all = [];
resp_all = [];
figure, hold on,
for subj = 1:10
    vec_acc = reshape(mean_acc_sel_len(:,:,subj),1,[]);
    vec_resp = reshape(resp.average_resp(:,:,subj),1,[]);
    if IS_ZSCORE
        vec_acc = zscore(vec_acc);
        vec_resp = zscore(vec_resp);
    end
    if IS_LINEFIT
        acc_all = [acc_all,vec_acc];
        resp_all = [resp_all,vec_resp];
    end
    for it = 1:length(vec_acc)
        scatter(vec_resp(it), vec_acc(it),50,cont_color{cont_list(it)},stim_symbol{stim_list(it)});
    end
end
if IS_LINEFIT
    [R,P] = corrcoef(resp_all,acc_all);
    [p,s,mu] = polyfit(resp_all,acc_all,1);
    yfit = p(1)*resp_all + p(2);
    plot(resp_all,yfit,'k-','linewidth',1)
    legend_tex = sprintf('y=%.2fx+%.2f',p(1),p(2));
    legend('FMSEQ', 'MSEQ', 'LF-SSVEP', 'HF-SSVEP',legend_tex,'location','northeastoutside');
else
    legend('FMSEQ', 'MSEQ', 'LF-SSVEP', 'HF-SSVEP', 'location','northeastoutside');
end
set(gcf,'position',[50, 50, 750, 500]); set(gca,'FontSize',16)

if IS_ZSCORE
    xlabel('Z-scored Ratings of Perceived Intensity')
    ylabel('Z-scored Decoding Accuracy')
    ylim([-2.5 2.5]); xlim([-2.5 2.5])
    if IS_LINEFIT
        export_fig('Fig5_scatter_zscore_fit','-png','-transparent');
    else
        export_fig('Fig5_scatter_zscore,'-png','-transparent');
    end
else
    xlabel('Ratings of Perceived Intensity')
    ylabel('Decoding Accuracy (%)')
    ylim([0 100]); xlim([1 7])
    if IS_LINEFIT
        export_fig('Fig5_scatter_raw_fit','-png','-transparent');
    else
        export_fig('Fig5_scatter_raw,'-png','-transparent');
    end
end

%% Fig 5-2: linear fitting curve to raw or z-scored data
IS_ZSCORE = 0;
acc_stim = zeros(4,30);
resp_stim = zeros(4,30);
for subj = 1:10
    vec_acc = reshape(mean_acc_sel_len(:,:,subj),1,[]);
    vec_resp = reshape(resp.average_resp(:,:,subj),1,[]);
    if IS_ZSCORE
        vec_acc = zscore(vec_acc);
        vec_resp = zscore(vec_resp);
    end
    for STIM = 1:4
        stim_idx = find(stim_list==STIM);
        acc_stim(STIM,(subj-1)*3+[1:3]) = vec_acc(stim_idx);
        resp_stim(STIM,(subj-1)*3+[1:3]) = vec_resp(stim_idx);
    end
end

legend_tex = cell(1,4);
figure,
for STIM = 1:4
    p = polyfit(resp_stim(STIM,:),acc_stim(STIM,:),1);
    yfit = p(1)*resp_stim(STIM,:) + p(2);
    scatter(resp_stim(STIM,:),acc_stim(STIM,:),15,cont_color{STIM},'filled'); hold on
    plot(resp_stim(STIM,:),yfit,'color',cont_color{STIM},'linewidth',1)
    legend_tex{STIM} = sprintf('y=%.2fx+%.2f',p(1),p(2));
end
legend('FMSEQ',legend_tex{1},'MSEQ',legend_tex{2},...
    'LF-SSVEP',legend_tex{3},'HF-SSVEP',legend_tex{4},'location','northeastoutside');
set(gcf,'position',[50, 50, 750, 500]); set(gca,'FontSize',16)
if IS_ZSCORE
    xlabel('Z-scored Ratings of Perceived Intensity')
    ylabel('Z-scored Decoding Accuracy')
    ylim([-2.5 2.5]); xlim([-2.5 2.5])
    export_fig('Fig5_linearfit','-png','-transparent');
else
    xlabel('Ratings of Perceived Intensity')
    ylabel('Decoding Accuracy (%)')
    ylim([0 100]); xlim([1 7])
    export_fig('Fig5_linearfit_raw','-png','-transparent');
end

