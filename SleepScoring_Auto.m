function OUT = SleepScoring_Auto(Filebase, params, CluParams, Path4Results)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% version information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 28 June 2017
%
% optimized parameters and post-adjustment procedures
%   1. smoothing iteration
%   2. applying new AW-REM violation rule
%
%% 26 June 2017
% fully automatic scoring by choosing REM and NREM clusters
%
%% 21 June 2017
% 1. added sorting by EMG
% 2. replaced EGM histograms with EMG box plot
%
%% 20 June 2017
% (almost) automatic sleep scoring

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HOW TO USE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filebase ... such as 'exp1' (BUT ASSUME 1kS .dat file, not .eeg or .dat@20kS)
% params
%   .ChInfo
%       .nch ... the total number of channels (e.g., 17)
%       .EEG ... EEG ch (e.g., 6)
%       .EMG ... EMG chs (single or bipolar) (e.g., 8 or [8 9])
%   .K ... # of clusters (default 5)
%   .Twin ... window size (default 4 sec)
%   .cWins ... the number of windows to smooth (default 4)

% CluParams
%   .Fs = 1000
%   .fpass = [0 45];
%   .tapers = [3 5];
%   .trialave = 0;
%   .max_iter = 500;
%   .lowcost = 2.5;
%   .theta = [6 9];
%   .EMGlim = 5;

% OUT
%   .Score (1, NREM; 2, REM; 3, AW)
%   .clus = Klus;
%   .pwr = Pwr; (freq x wins) in dB
%   .freq = f;
%   .emg = EMGrms;
%   .tsne = mappedX;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART I (fully automated process)
% 1-1. taking EEG signals and compute powerspectrogram (no overlapping window) (0-45 Hz)
% 1-2. dimensionality reduction with t-SNE
% 1-3. clustering with a Gaussian mixture model in 2 t-SNE dimensions and EMG rms (K = 5)
%   preclustering ... K-1 means clustering + putative REM assignment
% 1-4. clusters are sorted by EMG RMS mean
% 
% PART II (semi-automatic process)
% % 2-1. evaluate the clustering result and identify REM and NREM clusters
% % (sometime requires re-clustering to separate a REM cluster)

% new algorith (26 June 2017)
% 2-1-1. REM cluster? (peak of mean PSD at 6-9 Hz?)
% 2-1-2. NREM is the EMG minimum cluster except REM clu 

% 2-2. smoothing transition (when state shifts, the next cWins windows should
% be constant)
% 2-3. defragmentation (keep defragmented states constant) ... cWins
% 2-4. error correction (if AW->REM, AW becomes REM)
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% for future improvement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1-3. clustering ... EM algorithm can be implemented to have a better result
% 2-1. may be possible to take an automatic process by accessing the mean
% PSD profiles (should be easy to pick a REM cluster at least)
% 2-2 & 2-3 ... need to assess exactly how this improves the score

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% empirical experience to date
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. if there are significant chunk of REM sleep epochs, this provides a
% reliable result (easy to identify a REM cluster). Performance should be
% >85% with a good skill to identify REM and NREM clusters
% 2. Single channel EMG is fine.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% key parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Twin ... window size to assess (default = 4 sec)
% 2. K ... # of clusters (default = 8)
% 3. cWins ... # of consecutive windows to check smooth transition (default = 4; cWins x Twin sec)
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% requirements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. tSNE toolbox with some edits
% 2. 


close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(params, 'ChInfo') == 0
    error('ERROR: Specify params.ChInfo!!\n'); 
elseif isfield(params.ChInfo, 'nch') == 0 || ...
        isfield(params.ChInfo, 'EEG') == 0 || isfield(params.ChInfo, 'EMG') == 0
    error('ERROR: You need .ChInfo.nch, .EEG, and .EMG!! \n');
end
if isfield(params, 'K') == 0
    params.K = 5;
end
if isfield(params, 'Twin') == 0
    params.Twin = 4;
end
if isfield(params, 'cWins') == 0
    params.cWins = 4;
end

if nargin < 3
    CluParams.Fs = 1000;
    CluParams.fpass = [1 45];
    CluParams.tapers = [3 5];
    CluParams.trialave = 0;
    CluParams.max_iter = 500;
    CluParams.lowcost = 2.5;
    CluParams.theta = [6 9];
    CluParams.EMGlim = 5;
end

if nargin < 4
    Path4Results = cd;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART I ... fully automatic process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if fopen([Path4Results, '\', Filebase, '_sleepscoring_auto.mat']) == -1
    OUT = SleepScoring_Auto_Cluster(Filebase, params.ChInfo, params.K, params.Twin, CluParams);
    % OUT.clus = Klus;
    % OUT.pwr = Pwr;
    % OUT.freq = f;
    % OUT.emg = EMGrms;
    % OUT.tsne = mappedX;
else 
    load([Path4Results, '\', Filebase, '_sleepscoring_auto.mat']);
end
fclose('all');

% params.cWins = 3;

%%%%%%%%
%% reclustering
% OUT = MyReClustering(OUT, params, CluParams);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% auto REM & NREM clusters pick up and post-adjustments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REM clu?
CandPwr = mean(OUT.pwr(:, OUT.clus == 1), 2);
[~, MaxIdx] = max(CandPwr, [], 1);
MaxFreq = OUT.freq(MaxIdx);
idx = find(MaxFreq > CluParams.theta(1) & MaxFreq < CluParams.theta(2)); % theta peak 
if ~isempty(idx) % there is a REM clu
    Cinfo.REM = 1;
    Cinfo.NREM = MySearchNREM(OUT, 2:3);
else % REM rescue
    [OUT, Cinfo.REM] = MyREMrescue(OUT, params, CluParams);
    % Cinfo.REM = [];
    Cinfo.NREM = MySearchNREM(OUT, 1:3);
end

%% smoothing & defragmentation
Score = MySmoothingDefrag(OUT, Cinfo, params.cWins);

%% final figure
CO = 'bgr';
Lbl = {'NREM','REM','AW'};
figure(5);
for k = 1:3
    idx = find(Score == k);
    if ~isempty(idx)
        MyPwr = OUT.pwr(:,idx);

        subplot(2,4,k); %% spectrum
        plot(OUT.freq, MyPwr, [CO(k), '-']);hold on;
        plot(OUT.freq, mean(MyPwr,2), 'k-'); hold off;
        xlim(CluParams.fpass);
        ylim([-75 -15]);
        title(Lbl{k});

        subplot(2,4,k+4); %% image
        imagesc(1:length(idx),OUT.freq,MyPwr, [-65 -30]);
    end
end
%% EMG rms
subplot(244); % EMG rms
boxplot(OUT.emg, Score, 'notch', 'on');
xlabel('clu#');ylabel('EMG rms');
ylim([0 0.1]);
box off;
set(gca,'XTick',1:3,'XTickLabel',Lbl);
title('EMGrms');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data store
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.Cinfo = Cinfo;
OUT.Score = Score;
Log = [];
save([Path4Results, '\', Filebase, '_sleepscoring_auto.mat'],'OUT','params','CluParams', 'Cinfo');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% this section is for data with manual scores
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% check a pre-existing file and comparisons
% if exist('ManualSleepScore/sleepscoring_results.mat') ~= 0
%     Name = 'ManualSleepScore/sleepscoring_results.mat';
%     tmp = load(Name,'Score');
%     [HitRate, HitIdx, HitDetail] = MyLocalComp(Score, tmp.Score);
%     ManuScore = tmp.Score;
% elseif exist('manual sleep scoring/sleepscoring_results.mat') ~= 0
%     Name = 'manual sleep scoring/sleepscoring_results.mat';
%     tmp = load(Name,'Score');
%     [HitRate, HitIdx, HitDetail] = MyLocalComp(Score, tmp.Score);
%     ManuScore = tmp.Score;
% else % no pre-existing file
%     save('sleepscoring_results.mat', 'Score', 'Log');
%     HitRate = [];
%     HitIdx = [];
%     HitDetail = [];
%     ManuScore = [];
% end
% 
% OUT.HitRate = HitRate;
% OUT.HitIdx = HitIdx;
% OUT.HitDetail = HitDetail;
% OUT.ManuScore = ManuScore;
% fprintf([Filebase, ': ', num2str(HitRate), ' matched\n']);

%% DONE !!


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sub-functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% core function for smoothing and defrag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Score = MySmoothingDefrag(INPUT, Cinfo, cWins)

%% assumption
% there must be a NREM clu.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clus = INPUT.clus;

%% where is REM?
REMidx = [];
if ~isempty(Cinfo.REM)
    for c = 1:length(Cinfo.REM)
        REMidx = [REMidx; find(clus == Cinfo.REM(c))];
    end
end

%% where is NREM?
NREMidx = [];
for c = 1:length(Cinfo.NREM)
    NREMidx = [NREMidx; find(clus == Cinfo.NREM(c))];
end

%% initialize sleep score
Score = 3*ones(length(clus),1); % default is AW = 3
Score(REMidx) = 2;
Score(NREMidx) = 1;

Score0 = Score;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% before adjustment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4);subplot(511); %% spectrogram
imagesc(1:length(Score), INPUT.freq, INPUT.pwr, [-65 -30]);
ylabel('Hz'); title('EEG spectrogram');

subplot(512); %% just after clustering
bar(Score, 'k');
box off;
axis([0 length(Score)+1 0.5 3.5]);
ylabel('states');
set(gca,'YTick',1:3,'YTickLabel',{'NREM','REM','AW'});
title('just after clustering');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% update (core)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 1. smoothing transition (iterative process until violation disappears
Score1 = Score;
% i = 1;
while 1    
    idx = find(diff(Score(1:end-cWins+1)) ~= 0) + 1; % state changes
    for w = 1:length(idx)
        MySs = Score(idx(w):idx(w)+cWins-1);
        if sum(diff(MySs)) ~= 0 % do not accept the change
            Score(idx(w)) = Score(idx(w)-1);
        end
    end
   
    
    %% check
    [HitRate, ~, ~] = MyLocalComp(Score, Score1);
    Score1 = Score;
    if HitRate == 1
        break;
    end

end

[HitRate, HitIdx, ~] = MyLocalComp(Score, Score0);

figure(4);
subplot(513); %% just after smoothing
bar(Score,'r');hold on;
bar(HitIdx, Score(HitIdx), 'k'); hold off;
box off;
axis([0 length(Score)+1 0.5 3.5]);
ylabel('states');
set(gca,'YTick',1:3,'YTickLabel',{'NREM','REM','AW'});
title(['after smoothing, changed by ' , num2str(100*(1-HitRate)), ' %']);

%% 2. error correction
while 1
    [Score, nErrors] = MyErrorCorrection(Score);
    if nErrors == 0
        break;
    end
end

[HitRate, HitIdx, ~] = MyLocalComp(Score, Score0);

%% after all
subplot(514); %% just after defrag
bar(Score,'r');hold on;
bar(HitIdx, Score(HitIdx), 'k'); hold off;
box off;
axis([0 length(Score)+1 0.5 3.5]);
xlabel('win#'); ylabel('states');
set(gca,'YTick',1:3,'YTickLabel',{'NREM','REM','AW'});
title(['after all adjustments, changed by ' , num2str(100*(1-HitRate)), ' %']);

subplot(515); % emg
plot(INPUT.emg,'k');
xlim([1, length(INPUT.emg)]);
ylim([0 1.1*max(INPUT.emg)]);
box off;
title('EMG pwr');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sub-sub-functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [HitRate, HitIdx, HitDetail] = MyLocalComp(New, Org)

if length(Org) > length(New)
    Org = Org(1:length(New));
elseif length(Org) < length(New)
    New = New(1:length(Org));
end

tmp = Org - New;
HitIdx = find(tmp == 0);
HitRate = length(HitIdx)/length(New);

HitDetail = zeros(1,6); % hit rates & epoches in Org
for s = 1:3
    OrgIdx = find(Org == s);
    NewIdx = find(New == s);
    if ~isempty(NewIdx)
        HitDetail(s) = length(intersect(OrgIdx,NewIdx))/length(NewIdx);
    elseif ~isempty(OrgIdx)
        HitDetail(s) = 0;
    else
        HitDetail(s) = 1;
    end
    
    HitDetail(s+3) = length(OrgIdx);
end
        

%%
function NewScore = MyDefrag(Score, State, cWins)

Idx = find(Score == State);
D = diff(Idx);
idx = find(D <= cWins-1 & D > 1);
if ~isempty(idx)
    for i = 1:length(idx)
        Score(Idx(idx(i)):Idx(idx(i)+1)) = State;
    end
end
NewScore = Score;

%%
function [NewScore, nErrors] = MyErrorCorrection(Score)

%% find AW(3) to REM(2) transition
D = diff(Score);
idx1 = find(D == -1);
idx2 = find(Score == 3);
idx = intersect(idx1, idx2);
if ~isempty(idx)
    for i = 1:length(idx)
        if ~isempty(find(Score(max([idx(i)-10, 1]):idx(i))) == 1) % no NREM 10 wins before?
            Score(idx(i)) = 2;
        end
    end
end
nErrors = length(idx);
NewScore = Score;

%%%%%%%%%%%%%%%%%%
function NREMidx = MySearchNREM(OUT, Kidx)

Frange = [0 10];
Mpwr = zeros(length(Kidx),1);

for k = 1:length(Kidx)
    MeanPwr = mean(OUT.pwr(:, OUT.clus == Kidx(k)), 2);
    idx = find(OUT.freq >= Frange(1) & OUT.freq <= Frange(2));
    Mpwr(k) = sum(MeanPwr(idx));
end

[~, idx] = max(Mpwr);
NREMidx = Kidx(idx);

%%%%%%%%%%%%%%%%%%%%%%%%
function [OUT, REMidx] = MyREMrescue(OUT, params, CluParams)

% OUT.clus = Klus;
% OUT.pwr = Pwr; (freq x wins)
% OUT.freq = f;
% OUT.emg = EMGrms;
% OUT.tsne = mappedX;

Pwr = OUT.pwr;
Z = OUT.emg;
f = OUT.freq;

[~, MaxIdx] = max(Pwr, [], 1);
MaxFreq = f(MaxIdx);
idx1 = find(MaxFreq > CluParams.theta(1) & MaxFreq < CluParams.theta(2)); % theta peak 
% idx2 = find(Z < prctile(Z, CluParams.EMGlim)); % 5% percentile
idx2 = find(Z < prctile(Z, 4)); % x% percentile
if ~isempty(idx1)
    idx = intersect(idx1, idx2);
    OUT.clus(idx) = params.K + 1; % extra cluster for REM
    REMidx = params.K + 1; 
else
    REMidx = [];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OUT = MyReClustering(OUT, params, CluParams)

K = params.K;
mappedX = OUT.tsne;
EMGrms = OUT.emg;
Pwr = OUT.pwr;
f = OUT.freq;

%%%%
CO = jet(K);
X = mappedX(:,1); Y = mappedX(:,2); Z = log10(EMGrms');
Data = [(X-mean(X))./std(X), (Y-mean(Y))./std(Y), (Z-mean(Z))./std(Z)];

%% initial condition
% K-1 means clustering
Int = kmeans(Data, K-1);

% addition putative REM
[~, MaxIdx] = max(Pwr, [], 1);
MaxFreq = f(MaxIdx);
idx1 = find(MaxFreq > CluParams.theta(1) & MaxFreq < CluParams.theta(2)); % theta peak 
idx2 = find(Z < prctile(Z, CluParams.EMGlim)); % 5% percentile
figure(4);
if ~isempty(idx1)
    REMidx = intersect(idx1, idx2);

    % update initial condition
    Int(REMidx) = K;
    for k = 1:K
        idx = find(Int == k);
        MyPwr = Pwr(:,idx);

        figure(4);subplot(3,ceil((K+1)/3),k); %% image
        imagesc(1:length(idx),f,MyPwr, [-65 -30]);
        title(['cluster#',num2str(k)]);
    end
    drawnow
end

%% for EM
options = statset('MaxIter', 10000);
gmfit = fitgmdist(Data, K, 'CovarianceType', 'full', ...
     'SharedCovariance', false, 'Options', options, 'Start', Int);
Klus = cluster(gmfit, Data);

%% sort by EMG rms
MeanEMG = grpstats(EMGrms, Klus, {'mean'});
[~, Idx] = sort(MeanEMG);

NewKlus = Klus;
for k = 1:K
    MyK = Idx(k);
    Klus(NewKlus == MyK) = k;
end

OUT.clus = Klus;

