% ------------------------------------------------------------------------
%               ANOVA test for reponses
% ------------------------------------------------------------------------
contlist = [2 8 16];
sub_num = 10;
sum = zeros(4,3,10); % To store average response (among 40 trials for each condition) for each of the 10 subjects

for i = 1:sub_num       % Loop through all 10 subjects
    name = sprintf('respMat_s%d.mat',i);
    currentStruct = load(name);         % Load reesponse, but the cell is loaded as a struct
    currentMat = currentStruct.respMat; % Convert the loaded struct back to cell
    statMat = ones(4,3,40);             % To store all responses for one subject
    index_mat = ones(4,3);              % A index matrix indicating which trial (out of 40) is being inputed next
    for k = 1:520                               % Loop through all trials
        rating = currentMat{2}.rate_percept(k);
        mode = currentMat{2}.stimuli(k);
        cont = currentMat{2}.contrast(k);
        if mode == 5             % Separate the case of image
        else                     % The cases of all othe stimulis
            contindex = find(cont == contlist);
            statMat(mode,contindex,index_mat(mode,contindex)) = rating;
            index_mat(mode,contindex) = index_mat(mode,contindex) + 1;
        end
    end
    
    for n = 1:4
        for m = 1:3
            sum(n,m,i) = mean(statMat(n,m,:)); % Take mean of 40 trials for each condition and move to the "sum" matrix
        end
    end
end

row_len = 4; % Different colomns are different stimulis, which makes length of each row 4
col_len = length(contlist); % Different colomns are different contrasts, which makes length of each column 3xsubject_number
mat_2way = zeros(col_len*sub_num,row_len);

for row = 1:col_len
    for col = 1:row_len
        for sub_n = 1:sub_num
            row_num = 10*(row - 1)+sub_n;
            mat_2way(row_num,col) = sum(col,row,sub_n);
        end
    end
end

anova2(mat_2way,10);
%%
% ------------------------------------------------------------------------
%               paired t-test for reponses
% ------------------------------------------------------------------------
result_Cell = cell(4,1);
main_Data = zeros(10,1);
compare_Data = zeros(10,1);
for STIM = 1:4
    result_Mat = zeros(3,3);
    for CONT = 1:3
        count = 1;
        main_Data(:) = sum(STIM,CONT,:);
        for STIM_COMP = 1:4
            if STIM_COMP ~= STIM
                compare_Data(:) = sum(STIM_COMP, CONT, :);
                [~,result_Mat(count,CONT)] = ttest(main_Data,compare_Data);
                count = count+1;
            end
        end
    end
    result_Cell{STIM} = result_Mat;
end
%%
% ------------------------------------------------------------------------
%               ANOVA test for accuracies
% ------------------------------------------------------------------------
mean_acc_var = cell(4,3,10);
mu_ci_var = cell(4,3,10);
for SUB = 1:10
    name = sprintf('NSFVEP_s%d.mat', SUB);
    temp = load(name); 
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

[~,~,stats] = anova2(stat_Mat,10);
%%
% ------------------------------------------------------------------------
%               paired t-test for accuracies
% ------------------------------------------------------------------------
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
