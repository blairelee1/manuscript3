%Registration
%example function: registerRLS('fileName',5000,1,1,'cluster');
%threshold: check the pixel values inside the myograph wire of the roi
%image which should be around 20.

function registerRLS(fName,framesPerBlock,toThreshold,invertPolarity,procType)
%reads dimensions of the raw data
[data,sampling,timeStamps,sizeT]=readRLS(fName,0,1,[],'frame','uint8');
%saves the first image for roi selection
refImg=data;
[opt,met]=imregconfig('monomodal');

%writes a new file for registered rls file
pathWrite=[erase(fName,'.rls'),'_reg.rls'];
fileWriteId = fopen(pathWrite,'w');
for k=1:(30*1024)
    fwrite(fileWriteId,0 , 'uint64');
end
fseek(fileWriteId,0*1024,-1 );
fwrite(fileWriteId,uint64(size(data,2)), 'uint64');
fwrite(fileWriteId,uint64(size(data,1)), 'uint64');
fwrite(fileWriteId,uint64(sizeT), 'uint64');
fwrite(fileWriteId,uint64(sampling), 'uint64');
fseek(fileWriteId,30*1024,-1 );
fwrite(fileWriteId,timeStamps, 'uint64');
fwrite(fileWriteId,refImg,'uint8');

%the roi rectangle should include both arms of the wire 
figure(1)
imagesc(refImg);
title('Select ROI with a registration marker')
h=imrect;
pos = wait(h);
pos=round(pos);
x1=pos(1);
x2=pos(1)+pos(3);
y1=pos(2);
y2=pos(2)+pos(4);
refImg=refImg(y1:y2,x1:x2);

contrastMetrics=zeros(ceil(sizeT/framesPerBlock),2);

if toThreshold==1
    figure(1)
    imagesc(refImg)
    threshold=input('Threshold: ');
    refImg(refImg<threshold)=0;
    refImg(refImg>=threshold)=1;
    refImg=uint8(refImg);
end
if invertPolarity==1
    refImg=max(refImg(refImg(:)~=inf))-refImg;
end

tic
for blockIdx=1:1:ceil(sizeT/framesPerBlock)
    framesToRead=min(framesPerBlock,sizeT-((blockIdx-1)*framesPerBlock+1));
    
    [data,~,timeStamps,~]=readRLS(fName,(blockIdx-1)*framesPerBlock+1,framesToRead,[],'frame','uint8');
    subdata=data(y1:y2,x1:x2,:);
    metricImgBefore=mean(single(subdata),3);
    
    if toThreshold==1
        subdata(subdata<threshold)=0;
        subdata(subdata>=threshold)=1;
        subdata=uint8(subdata);
    end
    if invertPolarity==1
        subdata=max(subdata(subdata(:)~=inf))-subdata;
    end
    
    switch procType
        case 'cpu'
            for i=1:1:framesToRead
                img=squeeze(subdata(:,:,i));
                tform=imregtform(img,refImg,'translation',opt,met);
                data(:,:,i)=imwarpSame(data(:,:,i),tform);
                
%                 %debug code
%                 figure(1)
%                 subplot(1,2,1)
%                 imshowpair(refImg,img)
%                 img=imwarpSame(img,tform);
%                 subplot(1,2,2)
%                 imshowpair(refImg,img)
%                 drawnow
            end
        case 'cluster'
            parfor i=1:framesToRead
                img=squeeze(subdata(:,:,i));
                tform=imregtform(img,refImg,'translation',opt,met);
                data(:,:,i)=imwarpSame(data(:,:,i),tform);
            end
    end
    metricImgAfter=mean(single(data(y1:y2,x1:x2,:)),3);
    for i=1:1:framesToRead
        fwrite(fileWriteId,timeStamps(i), 'uint64');
        fwrite(fileWriteId,squeeze(data(:,:,i)),'uint8');
    end
    
    figure(2)
    subplot(1,2,1)
    imagesc(metricImgBefore)
    subplot(1,2,2)
    imagesc(metricImgAfter)
    title(['Block ',num2str(blockIdx)]);
    drawnow
    
    contrastMetrics(blockIdx,1)=std(metricImgBefore(:))./mean(metricImgBefore(:));
    contrastMetrics(blockIdx,2)=std(metricImgAfter(:))./mean(metricImgAfter(:));
    
    disp([num2str(blockIdx),' out of ',num2str(ceil(sizeT/framesPerBlock)),...
        ' blocks registered in ',num2str(toc),'s. Contrast metrics: ',...
        num2str(std(metricImgBefore(:))./mean(metricImgBefore(:))),' -before, ',...
        num2str(std(metricImgAfter(:))./mean(metricImgAfter(:))),' -after']);
end
fclose all;

figure
plot(contrastMetrics(:,1));
hold on
plot(contrastMetrics(:,2));
hold off
end


