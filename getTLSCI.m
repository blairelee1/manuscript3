%getTLSCI - calculates temporal Laser Speckle Contrast Images
%
% Syntax:  output1 = function_name(input1,input2,input3,input4)
%
% Inputs:
%    data       - raw laser speckle data as 3d [y,x,t] matrix
%    kernelSize - number of pixels in a side of the kernel
%    procType   - choose the processor type: use 'cpu' or 'gpu'
%    dsType     - downsampling type result is either same size as data or
%                 downsampled by kernel size. Use: 'none' or 'kernel'
%
% Outputs:
%    tLSCI      - processed data as [y,x,t] 3d matrix
%
% Example:
%    tLSCI=getSLSCI(data,25,'gpu','none')
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: getSLSCI.m

% Author: DD Postnov, PhD
% BOAS lab, Boston University
% email address: dpostnov@bu.edu
% Last revision: 3-March-2018

%------------- BEGIN CODE --------------

function tLSCI=getTLSCI(data,kernelSize,procType,dsType)
if strcmp(dsType,'none')
    T=1:1:size(data,3)-kernelSize+1;
elseif strcmp(dsType,'kernel')
    T=1:kernelSize:size(data,3)-kernelSize+1;
else
    T=1:1:size(data,3)-kernelSize+1;
end

tLSCI=zeros(size(data,1),size(data,2),length(T),'single');
counter=1;
if strcmp(procType,'cpu')
    for i=T
        frames=single(data(:,:,i:i+kernelSize-1));
        frameMean=squeeze(mean(frames,3));
        frameSTD=squeeze(std(frames,0,3));
        tLSCI(:,:,counter)=frameSTD./frameMean;
        counter=counter+1;
    end
elseif strcmp(procType,'gpu')
    for i=T
        frames=gpuArray(single(data(:,:,i:i+kernelSize-1)));
        frameMean=squeeze(mean(frames,3));
        frameSTD=squeeze(std(frames,0,3));
        tLSCI(:,:,counter)=gather(frameSTD./frameMean);
        counter=counter+1;
    end
    
elseif strcmp(procType,'fastcpu')
    frameSTD=movstd(single(data),[0,kernelSize-1],0,3);
    frameMean=movmean(single(data),[0,kernelSize-1],3);
    tLSCI=frameSTD./frameMean;
    tLSCI=tLSCI(:,:,T);
elseif strcmp(procType,'fastgpu') && strcmp(dsType,'none')
    gpu=gpuDevice;
    memoryToUse=0.1;
    dataTypeSize=4;
    pixelsN=floor((memoryToUse.*gpu.AvailableMemory)./(size(data,3).*dataTypeSize))
    if size(data,1)<=pixelsN
        blockSize=floor(pixelsN./(size(data,1)))
        blocksN=ceil(size(data,2)./blockSize)
        for i=1:1:blocksN
            i
            tic
            idxR=min(i*blockSize,size(data,2));
            idxL=(i-1)*blockSize+1;
            frames=gpuArray(single(data(:,idxL:idxR,:)));
            frameSTD=movstd(frames,[0,kernelSize-1],0,3);
            frameMean=movmean(frames,[0,kernelSize-1],3);
            tLSCI(:,idxL:idxR,:)=gather(frameSTD(:,:,(floor(kernelSize/2)+1):(end-floor(kernelSize/2)))./frameMean(:,:,(floor(kernelSize/2)+1):(end-floor(kernelSize/2))));
            toc
        end
    else
        error('Not sufficient memory to perform fast gpu operation');
    end
    
end

end

%------------- END OF CODE --------------
% Comments: large input data can lead to the memory overflow, particularily
% when uint8 imput data is provided. This can be controlled by additional
% outer loop and/or by conversion sLSCI data to scaled integer or by
% allowing downsampling by kernel size



