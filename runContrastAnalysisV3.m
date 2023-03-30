%% Clean up (optional)
clc
close all
clear all
%% Set parameters
%Compuational parameters
memoryLimit=100000; % RAM limit - set it to approx 60-80% of availible free RAM memory
lsciType="tLSCI"; %tLSCI or sLSCI
lsciKernel=25; % 5 for sLSCI, 25 for tLSCI
if strcmp(lsciType,'tLSCI')
    memoryLimit=memoryLimit/2;
end

%Output parameters
samplingRequest=1000; %for temporal averaging and downsampling
tdsKernel=0; %for temporal averaging and downsampling
sdsKernel=0; %for spatial downsampling



% File names: post raw registered 
fileNames{1}='Norm1LowBase.rls';
fileNames{2}='Norm2LowBase.rls';
fileNames{3}='Norm3LowBase.rls';
fileNames{4}='Norm4LowBase.rls';
fileNames{5}='Norm5LowBase.rls';
fileNames{6}='Norm6LowBase.rls';
fileNames{7}='Norm7LowBase.rls';
fileNames{8}='Norm9LowBase.rls';
fileNames{9}='Norm10LowBase.rls';

fileNames{10}='Norm1LowPost.rls';
fileNames{11}='Norm2LowPost.rls';
fileNames{12}='Norm3LowPost.rls';
fileNames{13}='Norm4LowPost.rls';
fileNames{14}='Norm5LowPost.rls';
fileNames{15}='Norm6LowPost.rls';
fileNames{16}='Norm7LowPost.rls';
fileNames{17}='Norm9LowPost.rls';
fileNames{18}='Norm10LowPost.rls';
disp('Parameters set')

%% Registration (motion artifacts reduction)


%% Basic contrast analysis
for i=1:1:length(fileNames)
    disp(['Processing file ',num2str(i),' out of ',num2str(length(fileNames))]);
    [data,sampling,timeStamps,sizeT]=readRLS(fileNames{i},0,1,[],'frame','uint8');
    sampling=double(sampling);
    sizeT=double(sizeT);    
    
    framesToAverage=max(1,max(round(samplingRequest/sampling),tdsKernel));
    
    sdsKernel=max(1,sdsKernel);
    sampling=framesToAverage*sampling;
    sizeT=floor(sizeT/framesToAverage);
    
    X=(1+floor(sdsKernel/2)):sdsKernel:(size(data,1)-floor(sdsKernel/2));
    Y=(1+floor(sdsKernel/2)):sdsKernel:(size(data,2)-floor(sdsKernel/2));
    
    switch lsciType
        case 'sLSCI'
            dataLSCI=zeros(length(X),length(Y),sizeT,'single');
            time=zeros(1,sizeT,'uint64');
        case 'tLSCI'
            dataLSCI=zeros(length(X),length(Y),sizeT-lsciKernel+1,'single');
            time=zeros(1,sizeT-lsciKernel+1,'uint64');
        otherwise
            disp('Error: wrong LSCI processing type')
            break;                       
    end
    
    
    
    memPerFrameInt8=size(data,1)*size(data,2)./1024./1024; %mbytes
    memPerFrameSingle=size(data,1)*size(data,2)./1025/1024*4; %mbytes
    
    blockSize=floor((memoryLimit-length(dataLSCI(:))/1024/1024*4)./(1.2*(memPerFrameInt8+memPerFrameSingle))./framesToAverage);
    procFramesPerBlock=blockSize;
    blockSize=floor(procFramesPerBlock.*framesToAverage);
    blocksN=ceil(sizeT/procFramesPerBlock);
    
    if blockSize<framesToAverage
        disp('Error: Spatial downsampling: negative or insufficient number of frames per processing block - check the memory restriction')
        break;
    end
    
    counter=0;
    for ii=1:1:blocksN
        tic
        switch lsciType
            case 'sLSCI'
                [data,~,timeStamps,~]=readRLS(fileNames{i},(ii-1)*blockSize,min(blockSize,sizeT*framesToAverage-(ii-1)*blockSize),[],'frame','uint8');
                data=getSLSCI(data,lsciKernel,'cpu','none');
            case 'tLSCI'
                [data,~,timeStamps,~]=readRLS(fileNames{i},(ii-1)*blockSize,min(blockSize+lsciKernel-1,sizeT*framesToAverage-(ii-1)*blockSize),[],'frame','uint8');
                if size(data,3)>=lsciKernel
                    data=getTLSCI(data,lsciKernel,'fastcpu','none');
                    
                    if size(data,3)~=min(blockSize,sizeT*framesToAverage-(ii-1)*blockSize-lsciKernel+1)
                        disp('Error: unexpected number of processed frames')
                        break
                    end
                else
                    disp('Error: insufficient number of frames')
                    break
                end
        end
        if floor(size(data,3)/framesToAverage)>0
            for iii=1:1:floor(size(data,3)/framesToAverage)
                img=data(:,:,((iii-1)*framesToAverage+1):iii*framesToAverage);
                time((ii-1)*procFramesPerBlock+iii)=round(mean(timeStamps(((iii-1)*framesToAverage+1):iii*framesToAverage)));
                if framesToAverage>1
                    img=squeeze(mean(img,3));
                end
                if sdsKernel>1
                    img=conv2(img,ones(sdsKernel),'same');
                    img=img./(sdsKernel^2);
                end
                dataLSCI(:,:,(ii-1)*procFramesPerBlock+iii)=img(X,Y);
            end
        end
        timeElapsed=toc;
        disp(['Finished block ',num2str(ii),' out of ',num2str(blocksN),'. Time elapsed ',num2str(timeElapsed)]);
    end
    save([ erase(fileNames{i},'.rls'),'.mat'],'dataLSCI','time','sampling','-v7.3');
end

