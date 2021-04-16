% ------------------------------------------------------------------------
%               TRCA-based VEP detection algorithm
% -------------------------------------------------------------------------
function [mean_acc, mu_ci] = vep_trca_var(epochedEEG,STIM,CONT,DATA_LEN_LIST,CH_LIST)

fprintf('Results of the ensemble TRCA-based method.\n');

if isempty(DATA_LEN_LIST)
    DATA_LEN_LIST = 1.2; % use 1.2 second 
end
if isempty(CH_LIST)
    % Selection: {[12,19],[18:20],[12,18:20],[12,17:21,24:26],[5,9:15,17:21,24:26],[5,9:15,16:32],[1:32]}; 
    CH_LIST = [12,17:21,24:26]; % 9 channels centered at Oz
end
    
% define settings 
numTrial = 10;
SRATE = 512;
plot_fig = 0;

% define TRCA parameters
is_ensemble = 1;    % use ensemble classifier (one classifier for each code seq)
DATA_OFFSET = 0.2+0.1;   % 0.2 sec baselien + 0.1 sec offset (estimated delay)
% CH_label = [1,2,3,4,9,16,25,32];


% define condition for testing
if length(DATA_LEN_LIST) > 1
    var_list = DATA_LEN_LIST;
else
    var_list = CH_LIST;
end


% main loop - cross validation for each condition
mean_acc = zeros(1, length(var_list));
mu_ci = zeros(length(var_list),2);

for var_i = 1:length(var_list)
    
    if length(DATA_LEN_LIST) > 1
        data_len = var_list(var_i);
        CH = CH_LIST;
    else
        data_len = DATA_LEN_LIST;
        CH = var_list{var_i};
    end

    data_range = floor(DATA_OFFSET*SRATE) + (1:1:floor(data_len*SRATE));
    NSAMP = length(data_range);
    
    % extract EEG epochs
    trial_eeg = zeros(4,length(CH),NSAMP,numTrial);
    for loc_i = 1:4
        for tr_i = 1:10
            trial_eeg(loc_i,:,:,tr_i) = epochedEEG{STIM,CONT,loc_i,tr_i}(CH,data_range);
        end
    end
    
    % Leave-one-trial-out cross validation classification accuracy
    labels = 1:4;
    accs = zeros(1,numTrial);
    for cv_i = 1:1:numTrial
        
        % Training stage
        traindata = trial_eeg;
        traindata(:, :, :, cv_i) = [];
        model = train_trca(traindata);
        
        % Test stage
        testdata = squeeze(trial_eeg(:, :, :, cv_i));
        estimated = test_trca(testdata, model, is_ensemble);
        
        % Evaluation
        is_correct = (estimated==labels);
        accs(cv_i) = mean(is_correct)*100;
        
    end % loocv_i
    
    % Summarize results
    alpha_ci = 0.05;    % alpha for computing confidence interval
    ci = 1-alpha_ci;    % confidence interval
    [mu, ~, muci, ~] = normfit(accs, alpha_ci);
    fprintf('Mean accuracy = %2.2f %% (%2d%% CI: %2.2f - %2.2f %%)\n',...
        mu, ci, muci(1), muci(2));
    
    mean_acc(var_i) = mu;
    mu_ci(var_i,:) = muci;

end

% plot classification results vs. number of channels for each contrast
if plot_fig
    plot(CH_label, mean_acc,'linewidth',2);
    xlabel('Training data legnth (sec)'); ylabel('Cross validation accuracy (''%)');
    legend('FMC','MSEQ','LF_SSVEP','HF_SSVEP');
    set(gca,'fontsize',12)
    ylim([0 100]);
end

end
