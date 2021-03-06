%% addpath 
addpath(genpath('C:\Users\shawn\Desktop\NSF_BIC\expVEP_PTB\'));
eeglab


%% load EEG and marker streams from XDF file
% filename = 'data\expVEP_050119.xdf';
[filename, pathname] = uigetfile('*.xdf');

% load and integrate EEG and marker streams
EEG = pop_loadxdf([pathname, filename], 'streamtype', 'EEG', 'exclude_markerstreams', {});
EEG = pop_select( EEG,'nochannel',{'Trig1'});

% select EOG and EEG channels
EOG = pop_select( EEG,'channel',35:37); % EX3-EX5
EEG = pop_select( EEG,'channel',1:34);  % A1-A32, EX1-EX2
disp({EEG.chanlocs.labels});
pop_saveset(EEG,'filename','exp_vep_s02','filepath',pathname);


%% load other data streams from XDF file
[Streams,FileHeader] = load_xdf([pathname, filename]);
for it = 1:length(Streams), disp(Streams{it}.info.name); end

disp(Streams{2}.info)
srate = str2double(Streams{2}.info.nominal_srate);
figure, plot(Streams{2}.time_series(1,60*srate:(120*srate+20*srate)))

% S01: 91.0422 sec (rising edge of photo-diode) vs. 91.0703 (event marker of stim onset)
%      photo diode leads 28.1 msec than event marker

% S02: 80.8533 sec (rising edge of photo-diode) vs. 80.8965 (event marker of stim onset)
%      photo diode leads 43.2 msec than event marker


%% Eye blink analysis for event synchronization
EOG = pop_eegfiltnew(EOG, 0.5,40,3380,0,[],1);     EEG = eeg_checkset( EEG );
EOG = pop_epoch( EOG, {  'Blink'  }, [-0.8         0.8], 'epochinfo', 'yes');
EOG = eeg_checkset( EOG );
     
Y = zeros(3,EOG.trials);
I = zeros(3,EOG.trials);
for trial_id = 1:EOG.trials
    [Y(:,trial_id), I(:,trial_id)] = max(squeeze(EOG.data(:,:,trial_id)),[],2);
end
channel = 3;
figure, hist(EOG.times(I(channel,:)))
median(EOG.times(I),2)
std(EOG.times(I),[],2)

% S01: peak of blink activity in EEG happened 95 to 100 msec after event marker, though variance is large
% S02: peak of blink activity in EEG happened 25 to 32 msec before event marker, though variance is large
% Using eye blink for syncrhonization is flawed given variability in
% subjects' prediction and long duration of blink (0.5 sec)


%% ------------------------------------------------------------------------
%                       EEG Preprocessing
% -------------------------------------------------------------------------
filename = 'exp_vep_s02.set';
filepath = 'data\';
EEG = pop_loadset([filepath filename]);

% high pass filtering
lowcutoffTh = 1.0;  % Hz
EEG = pop_eegfiltnew(EEG,[],lowcutoffTh,[],1,[],0);

% % remove bad channels
% EEG = clean_flatlines(EEG);
% EEG = clean_channels_nolocs(EEG);
% 
% % remove line noise using cleanLine
% EEG = pop_cleanline(EEG, 'bandwidth',2,'chanlist', 1:EEG.nbchan,'computepower',1,'linefreqs',...
%     [60 120] ,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',1,'sigtype',...
%     'Channels','tau',100,'verb',0,'winsize',4,'winstep',1, 'ComputeSpectralPower','False');
% 
% % repair artifacts using ASR with mild threshold
% asr_stdcutoff = 20;
% EEG = clean_asr(EEG,asr_stdcutoff);

% store pre-processed EEG
EEG = eeg_checkset(EEG);
pop_saveset(EEG);


%% ------------------------------------------------------------------------
%               Extract EEG epochs and save data matrix
% -------------------------------------------------------------------------

% load preprocessed EEG data
filename = 'expVEP_050119_prep.set';
filepath = 'data\';
EEG = pop_loadset([filepath filename]);

% load behavioral results and stimulation codes
load('data\respMat_050119.mat')
code = respMat{1};

% extract EEG epoch
epochWindow = [-0.2, 3.2];   % define time window for EEG epoch (in sec)
numTrial = 10;
expression = 'Stim(?<stim>\d+)_CT(?<cont>\d+)_Loc(?<loc>\d+)';  % decode event markers
exp1_end_event_id = find(strcmp({EEG.event.type},'Exp2'))-1;




% % first experiment
% stimList = {'00','10','20','30','50'};   % stimulation pattern: 00: fmc-VEP (independent) 01: fmc-VEP (inverted) 02: fmc-VEP (chrome) 10: mseq
% contList = {'2','4','8'};          % contrast level: +- 4, 8, 64
% locList = {'1','2'};                % location: 1: left screen, 2: right screen
% 
% epochedEEG_exp1 = cell(length(stimList),length(contList),length(locList),numTrial);
% epochedIndex_exp1 = zeros(length(stimList),length(contList),length(locList));
% for event_it = 1:exp1_end_event_id
%     tokenNames = regexp(EEG.event(event_it).type,expression,'names');
%     if ~isempty(tokenNames)
%         stimIndex = find(strcmp(tokenNames.stim,stimList));
%         contIndex = find(strcmp(tokenNames.cont,contList));
%         locIndex = find(strcmp(tokenNames.loc,locList));
%         epochedIndex_exp1(stimIndex,contIndex,locIndex) = epochedIndex_exp1(stimIndex,contIndex,locIndex) + 1;
%         
%         dataIndex = floor(EEG.event(event_it).latency) + floor(epochWindow(1)*EEG.srate) : ...
%             floor(EEG.event(event_it).latency) + floor(epochWindow(2)*EEG.srate);
%         epochedEEG_exp1{stimIndex,contIndex,locIndex,epochedIndex_exp1(stimIndex,contIndex,locIndex)} = EEG.data(:,dataIndex);
%     end
% end
% save('data\S3_epochedEEG_exp1.mat','epochedEEG_exp1');

% % second experiment 
% locList = {'20','30','35','40','45','50','55','60','65','70','80'};                % location: 1: left screen, 2: right screen
% epochedEEG_exp2 = cell(length(locList),numTrial);
% epochedIndex_exp2 = zeros(1,length(locList));
% for event_it = exp1_end_event_id+1 : length(EEG.event)
%     tokenNames = regexp(EEG.event(event_it).type,expression,'names');
%     if ~isempty(tokenNames)
%         locIndex = find(strcmp(tokenNames.loc,locList));
%         epochedIndex_exp2(locIndex) = epochedIndex_exp2(locIndex) + 1;
%         
%         dataIndex = floor(EEG.event(event_it).latency) + floor(epochWindow(1)*EEG.srate) : ...
%             floor(EEG.event(event_it).latency) + floor(epochWindow(2)*EEG.srate);
%         epochedEEG_exp2{locIndex,epochedIndex_exp2(locIndex)} = EEG.data(:,dataIndex);
%     end
% end
% save('data\S3_epochedEEG_exp2.mat','epochedEEG_exp2');


%% Archive
% % apply high pass fileter
% highcutoffTh = 50;  % Hz
% EEG = pop_eegfiltnew(EEG,[],CONFIG.filter_lp_cutoff,[],0,[],plotfreqz);
% 
% 
% lpFilt = designfilt('lowpassiir', 'FilterOrder', 8, ...
%     'PassbandFrequency', highcutoffTh, 'PassbandRipple', 0.2,...
%     'SampleRate', EEG.srate);
% 
% tempData = zeros(EEG.nbchan,EEG.pnts);
% for ch_id = 1:EEG.nbchan 
%     tempData(ch_id,:) = filtfilt(lpFilt, tempData(ch_id,:));
%     fprintf('.');
% end
% 
% % visualize filter and filtered signals
% % fvtool(d3,'Fs',EEG.srate)
% 
% [praw,fraw] = periodogram(EEG.data(1,:),[],[],EEG.srate);
% [pbutt,fbutt] = periodogram(tempData(1,:),[],[],EEG.srate);
% figure, plot(fraw, 20*log10(abs(praw)),'b')
% hold on, plot(fbutt,20*log10(abs(pbutt)),'r')

% d1 = designfilt('bandstopiir','FilterOrder',2, ...
%     'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
%     'DesignMethod','butter','SampleRate',EEG.srate);
% d2 = designfilt('bandstopiir','FilterOrder',2, ...
%                 'HalfPowerFrequency1',119,'HalfPowerFrequency2',121, ...
%                 'DesignMethod','butter','SampleRate',EEG.srate);
% d3 = designfilt('bandstopiir','FilterOrder',2, ...
%                 'HalfPowerFrequency1',179,'HalfPowerFrequency2',181, ...
%                 'DesignMethod','butter','SampleRate',EEG.srate);
% d4 = designfilt('bandstopiir','FilterOrder',2, ...
%                 'HalfPowerFrequency1',239,'HalfPowerFrequency2',241, ...
%                 'DesignMethod','butter','SampleRate',EEG.srate);
%     tempData(ch_id,:) = filtfilt(d1, double(EEG.data(ch_id,:)));
%     tempData(ch_id,:) = filtfilt(d2, tempData(ch_id,:));
%     tempData(ch_id,:) = filtfilt(d3, tempData(ch_id,:));
%     tempData(ch_id,:) = filtfilt(d4, tempData(ch_id,:));


