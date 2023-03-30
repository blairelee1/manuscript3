%% SD file names 
clear 
clc
fNames{1} = 'SHR1HighBase_reg_segmented.mat';
fNames{2} = 'SHR1HighPost_reg_segmented.mat';
fNames{3} = 'SHR2HighBase_reg_segmented.mat';
fNames{4} = 'SHR2HighPost_reg_segmented.mat';
fNames{5} = 'SHR3HighBase_reg_segmented.mat';
fNames{6} = 'SHR3HighPost_reg_segmented.mat';
fNames{7} = 'SHR4HighBase_reg_segmented.mat';
fNames{8} = 'SHR4HighPost_reg_segmented.mat';
fNames{9} = 'SHR5HighBase_reg_segmented.mat';
fNames{10} = 'SHR5HighPost_reg_segmented.mat';
fNames{11} = 'SHR6HighBase_reg_segmented.mat';
fNames{12} = 'SHR6HighPost_reg_segmented.mat';
fNames{13} = 'SHR7HighBase_reg_segmented.mat';
fNames{14} = 'SHR7HighPost_reg_segmented.mat';
fNames{15} = 'SHR8HighBase_reg_segmented.mat';
fNames{16} = 'SHR8HighPost_reg_segmented.mat';
fNames{17} = 'SHR9HighBase_reg_segmented.mat';
fNames{18} = 'SHR9HighPost_reg_segmented.mat';


%% Analyize TGF signal and add parameters 

for i = 1:length(fNames)
i
load(fNames{i});

% ASLT parameters 
fois=0.01:0.0005:0.05;
Ncyc=15;
ord=[1,50];
fs=1;
% bandpass filter parameters 
fpass = [0.015 0.04]; % passband frequency 
% auc calculation parameters 
flimits = (fois >= fpass(1)) & (fois <= fpass(2));
% peak limit denominator 
denom = 2;

% add parameter to results 
result.fois = fois; 
result.Ncyc = Ncyc;
result.ord = ord;
result.fs = fs;
result.fpass = fpass;



% calculate mean bfi and aslt 
for ii = 1:length(result.segmves)

    %%%%%% prepare signal %%%%%%
    % tsClean : clean signal 
    tsClean = result.segmves(ii).K;
    tsClean = mean(tsClean,1,'omitnan');
    tsClean = 1./tsClean.^2;
    tsClean(isinf(tsClean)) = [];
    % tsPass: filtered signal 
    tsPass = bandpass(tsClean,fpass,fs);
    % tmp: calculate tmp on tsClean
    tmpNoisy = aslt(tsClean, fs, fois, Ncyc, ord, 0);
    tmpNoisy = squeeze(mean(tmpNoisy,2,'omitnan'));
    % tmpPass: calculate tmp on tsClean
    tmpPass = aslt(tsPass, fs, fois, Ncyc, ord, 0);
    tmpPass = squeeze(mean(tmpPass,2,'omitnan'));
    
    
    %%%%%%% collect parameters %%%%%%%
    meanBfi = mean(tsClean,'omitnan');
    sigma = std(tsPass,'omitnan'); 
    aucNoisy = trapz(tmpNoisy(flimits));
    aucPass = trapz(tmpPass(flimits));
    SNRauc = aucPass/(aucNoisy-aucPass);
    SNRsig = snr(tsClean,1); %
    [~, ~, ~, p] = findpeaks(tmpPass,fois); % find all the peaks 
    lim = std(p)/denom;
    [pks, locs, w, p] = findpeaks(tmpPass,fois,'MinPeakProminence',lim,'SortStr','descend'); % peaks with minimum prominence limit 
    numPks = length(pks);

    %%%%%% save analyzed signal and parameters to results %%%%%%%%%
    result.segmves(ii).tsClean = tsClean;
    result.segmves(ii).tsPass = tsPass;
    result.segmves(ii).tmpNoisy = tmpNoisy;
    result.segmves(ii).tmpPass = tmpPass;
    result.segmves(ii).meanBfi = meanBfi;
    result.segmves(ii).sigma = sigma;
    result.segmves(ii).aucNoisy = aucNoisy;
    result.segmves(ii).aucPass = aucPass;
    result.segmves(ii).SNRsig = SNRsig; % signal SNR expressed in decibals 
    result.segmves(ii).SNRauc = SNRauc; % auc SNR 
    result.segmves(ii).numPks = numPks;
    result.segmves(ii).pks = pks; % power 
    result.segmves(ii).locs = locs; % frequency 
    result.segmves(ii).w = w; % half prominence width 
    result.segmves(ii).p = p; % prominence
end 
save(fNames{i},'result','-v7.3')
end 




