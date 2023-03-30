%% register masks 
% register and get transform
% fixed image: post mask 
% "moving" image: baseline mask 
% the goal is to see how much the base mask (moving) has to move to match
% the mask location for post mask (fixed). 
%set file names
baseContrast='SHR8HighBase_reg.mat';
baseSegmented='SHR8HighBase_reg_segmented.mat';
postContrast='SHR8HighPost_reg.mat';
postSegmented='SHR8HighPost_reg_segmented.mat';


%% first check if registration is needed
load(baseSegmented)
baseResult=result;
baseBfiImg=baseResult.meanBfi;
baseMask=baseResult.maskId;
baseMask(baseMask ~= 0) = 1;
baseMask = logical(single(baseMask));
vesCount=length(baseResult.segmves);
% switch maskIni to maskId since that's the mask after manual segmentation 
load(postSegmented);
postResult=result;
postBfiImg=postResult.meanBfi;
postMask=postResult.maskId;
postMask(postMask ~= 0) = 1;
postMask = logical(single(postMask));
vesCount2=length(postResult.segmves);


figure
subplot(1,2,1)
imagesc(baseBfiImg)
caxis([prctile(baseBfiImg(:),5) prctile(baseBfiImg(:),99)]);
axis image
hold on 
visboundaries(baseMask)
hold off

subplot(1,2,2)
imagesc(postBfiImg)
caxis([prctile(postBfiImg(:),5) prctile(postBfiImg(:),99)]);
axis image
hold on 
visboundaries(postMask)
hold off
clearvars result
%% control point based registration 
moving=baseResult.maskId;
fixed=postResult.maskId;
[mp,fp] = cpselect(moving,fixed,Wait=true); % select control points interactively. Close the window once complete 

%% infer geometric transformation 
t = fitgeotform2d(mp,fp,'similarity');
Rfixed = imref2d(size(fixed));
registered = imwarp(moving,t,OutputView=Rfixed);
imshowpair(fixed,registered,"blend")

%% preview mask to post data to ensure alignment 
imagesc(postBfiImg)
caxis([prctile(postBfiImg(:),5) prctile(postBfiImg(:),99)]);
hold on 
visboundaries(registered)
hold off 
%%
% get config information from baseline segmented data
baseconfig=baseResult.config;
% file to apply the new mask 
config.fName=postContrast;
% set config parameters to baseline data 
config.toMaskMean=baseconfig.toMaskMean;
config.toMaskMin=baseconfig.toMaskMin;
config.maskErode=baseconfig.maskErode;

config.mfiltSize=baseconfig.mfiltSize;
config.gfiltSize=baseconfig.gfiltSize;
config.sfiltConf=baseconfig.sfiltConf;
config.atSensitivity=baseconfig.atSensitivity;
config.atNeighSize=baseconfig.atNeighSize;
config.prctiles=baseconfig.prctiles;

config.binarizeType=baseconfig.binarizeType;
config.relVal=baseconfig.relVal;

config.bwDilateErodeSize=baseconfig.bwDilateErodeSize;
config.bwareafiltSize=baseconfig.bwareafiltSize;
T=baseResult.aT;
clearvars result

%% load post registered contrast data 

load(config.fName)
config.sampling=sampling;
config.fS=1000./config.sampling;
time=(0:size(dataLSCI,3)-1)./config.fS;
img=mean(dataLSCI,3,'omitnan');
meanBfi=1./(img.^2);
imgToSegm=sqrt(meanBfi);
% figure
% img=meanBfi;
% imagesc(img)
% caxis([prctile(img(:),5),prctile(img(:),99)])
% title('BFI')
% imagesc(squeeze(mean(dataLSCI,3,'omitnan')));
%% final extraction using baseline registered mask  

maskIni=registered;
maskFinal=registered;

dataLSCI=reshape(dataLSCI,size(meanBfi,1)*size(meanBfi,2),length(time));
result.time=time;
result.config=config;
result.maskIni=maskIni;
result.meanBfi=meanBfi;
result.imgToSegm=imgToSegm;
result.aT=T;

rprops = regionprops(maskFinal,'PixelList','PixelIdxList','Centroid');
result.maskId=zeros(size(maskFinal));
for i=1:length(rprops)
    result.maskId(rprops(i).PixelIdxList)=i;
    result.segmves(i).rprops=rprops(i);
    result.segmves(i).id=i;
    result.segmves(i).rprops.area=length(rprops(i).PixelIdxList);
    result.segmves(i).K=dataLSCI(rprops(i).PixelIdxList,:);
    result.segmves(i).tsK=squeeze(mean(result.segmves(i).K,1));
%   result.segmves(i).tsBFI=1./(result.segmves(i).tsK.^2);
    ts = 1./(result.segmves(i).tsK.^2);
    ts(isinf(ts)) = []; % fix tail end data 
    result.segmves(i).tsBFI = ts;
    
    if sum(ismissing(ts)) > 0 % calculate ph and pow only if the vessel signal is entirely is preserved 
        result.segmves(i).ph=[];
        result.segmves(i).pow=[];
        result.segmves(i).frq=[];
    else
    [wt,f] = cwt(result.segmves(i).tsBFI,config.fS);
    result.segmves(i).ph=angle(wt);
    result.segmves(i).pow=single(abs(wt));
    result.segmves(i).frq=round(f,4);
    end 
end


save(strrep(config.fName,'.mat','_segmentedv2.mat'),'result','-v7.3');

