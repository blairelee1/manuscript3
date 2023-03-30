%% Settings
close all
clear all
%%
config.fName='SHR1HighBase_reg.mat';

config.toMaskMean=[0.1,0.20]; %check average contrast image and select values between myograph and high values.
config.toMaskMin=0.01;
config.maskErode=[51,3];

config.mfiltSize=[3,3];
config.gfiltSize=1;
config.sfiltConf=[3,1];
config.atSensitivity=0.48; %adjust depending on segmentation but typically .5 is sufficient 
config.atNeighSize=101;
config.prctiles=[0.1,99.9];

config.binarizeType='imbinarize';
config.relVal=0.995; % to be used together with config.binarizeType='config.relVal'

config.bwDilateErodeSize=3;
config.bwareafiltSize=[400,10000];
 


%% Load data
load(config.fName)
config.sampling=sampling;
config.fS=1000./config.sampling;
time=(0:size(dataLSCI,3)-1)./config.fS;
img=mean(dataLSCI,3,'omitnan');
meanBfi=1./(img.^2);
imgToSegm=sqrt(meanBfi);
figure
img=meanBfi;
imagesc(img)
caxis([prctile(img(:),5),prctile(img(:),99)])
title('BFI')

%% if sampling is 500 instead of 1000 use this section to resample data 
%with sampling that is 1000, skip this section!
if sampling ==500
dataLSCI(:,:,1:2:end-1)= (dataLSCI(:,:,1:2:end-1)+ dataLSCI(:,:,2:2:end))./2;
dataLSCI(:,:,2:2:end)=[];
sampling=1000;
end
%% Automated masking and segmentation
mask=dataLSCI>config.toMaskMean(1) & dataLSCI<config.toMaskMean(2);
mask=round(mean(mask,3));
mask=imerode(mask,ones(config.maskErode(1)));
%mask=mask & min(dataLSCI>config.toMaskMin,[],3) & min(~isnan(dataLSCI),[],3);
mask=imerode(mask,ones(config.maskErode(2)));
figure
imagesc(mask)
 
%%
extSize=floor(config.atNeighSize/2)+1;
fImgToSegm=imgToSegm;
fImgToSegm=medfilt2(fImgToSegm,config.mfiltSize);
fImgToSegm=imgaussfilt(fImgToSegm,config.gfiltSize);
valLow=double(prctile(fImgToSegm(mask(:)>0),config.prctiles(1)));
valHigh=double(prctile(fImgToSegm(mask(:)>0),config.prctiles(2)));
fImgToSegm=mat2gray(fImgToSegm,[valLow,valHigh]);
fImgToSegm=imsharpen(fImgToSegm,'Radius',config.sfiltConf(1),'Amount',config.sfiltConf(2));
fImgToSegm=mat2gray(fImgToSegm);
fImgToSegm=padarray(fImgToSegm,[extSize,extSize],'symmetric' );
T=adaptthresh(fImgToSegm,config.atSensitivity,'NeighborhoodSize',config.atNeighSize,'Statistic','median');
switch config.binarizeType
    case 'imbinarize'
        bw1=imbinarize(fImgToSegm,T);
    case 'config.relVal'
        bw1=(fImgToSegm./T)>config.relVal;
end
bw1=bw1(extSize+1:end-extSize,extSize+1:end-extSize);
fImgToSegm=fImgToSegm(extSize+1:end-extSize,extSize+1:end-extSize);
T=T(extSize+1:end-extSize,extSize+1:end-extSize);
bw1=bw1 & mask;
bw2 = imerode(bw1,ones(config.bwDilateErodeSize));
bw2 = imdilate(bw2,ones(config.bwDilateErodeSize));
maskIni=bwareafilt(bw2,config.bwareafiltSize,4);
maskFinal=maskIni;
figure
subplot(2,3,1)
img=meanBfi;
imagesc(img)
caxis([prctile(img(:),5),prctile(img(:),99)])
title('BFI')
subplot(2,3,2)
img=fImgToSegm;
imagesc(img)
caxis([prctile(img(:),5),prctile(img(:),99)])
title('Image to segment')
subplot(2,3,3)
img=fImgToSegm./T;
threshNormImg = img;
imagesc(img)
caxis([prctile(img(:),5),prctile(img(:),99)])
title('Threshold normalized')
subplot(2,3,4)
imagesc(bw1)
title('Thresholded mask')
subplot(2,3,5)
imagesc(maskIni)
title('Cleaned up mask')
subplot(2,3,6)
img=fImgToSegm./T;
imagesc(img)
caxis([prctile(img(:),5),prctile(img(:),99)])
hold on
visboundaries(maskIni,'Color','m');
hold off
title('Final result')
%% Manual mask correction - remove
close(figure(1));
h=figure(1);
img=fImgToSegm./T;
imagesc(img)
caxis([prctile(img(:),5),prctile(img(:),99)])
hold on
visboundaries(maskFinal,'Color','m');
hold off
xticks([]);
yticks([]);
sgtitle('Press spacebar to select a region to delete. Close the figure to exit.')
while true
    waitfor(h, 'CurrentCharacter', char(32))
    if exist('h') && ishandle(h)
        set(h, 'CurrentCharacter', char(12))
        correctionMask=roipoly;
        maskFinal(correctionMask)=0;        
        figure(1)
        img=fImgToSegm./T;
        imagesc(img)
        caxis([prctile(img(:),5),prctile(img(:),99)])
        hold on
        visboundaries(maskFinal,'Color','m');
        hold off
        xticks([]);
        yticks([]);        
    else
        break;
    end
end


%% Manual mask correction - add
close(figure(1));
h=figure(1);
img=fImgToSegm./T;
imagesc(img)
caxis([prctile(img(:),5),prctile(img(:),99)])
hold on
visboundaries(maskFinal,'Color','m');
hold off
xticks([]);
yticks([]);
sgtitle('Press spacebar to select a region to add. Close the figure to exit.')
while true
    waitfor(h, 'CurrentCharacter', char(32))
    if exist('h') && ishandle(h)
        set(h, 'CurrentCharacter', char(12))
        correctionMask=roipoly;
        maskFinal(correctionMask)=1;        
        figure(1)
        img=fImgToSegm./T;
        imagesc(img)
        caxis([prctile(img(:),5),prctile(img(:),99)])
        hold on
        visboundaries(maskFinal,'Color','m');
        hold off
        xticks([]);
        yticks([]);        
    else
        break;
    end
end


%% Final extraction of segmented vessels charachteristics
dataLSCI=reshape(dataLSCI,size(meanBfi,1)*size(meanBfi,2),length(time));
result.time=time;
result.config=config;
result.maskIni=maskIni;
result.meanBfi=meanBfi;
result.imgToSegm=imgToSegm;
result.aT=T;


rprops = regionprops(maskFinal,'PixelList','PixelIdxList','Centroid');
result.maskId=zeros(size(maskFinal));
for i=1:1:length(rprops)
    result.maskId(rprops(i).PixelIdxList)=i;
    result.segmves(i).rprops=rprops(i);
    result.segmves(i).id=i;
    result.segmves(i).rprops.area=length(rprops(i).PixelIdxList);
    result.segmves(i).K=dataLSCI(rprops(i).PixelIdxList,:);
    result.segmves(i).tsK=squeeze(mean(result.segmves(i).K,1));
    result.segmves(i).tsBFI=1./(result.segmves(i).tsK.^2);
    
end



save(strrep(config.fName,'.mat','_segmented.mat'),'result','-v7.3');