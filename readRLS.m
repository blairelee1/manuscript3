function [data,sampling,timeStamps,sizeT]=readRLS(fileName,startT,framesN,ROI,type,dataType)
C = strsplit(fileName,'.');
if ~strcmp(C{end},'rls')
   fileName=[fileName,'.rls']; 
end
counter=0;
fileReadId = fopen(fileName, 'r');
fseek(fileReadId,0*1024,-1 );
%NOTE IN SOME RLS FILES Y and X were misplaced!!! Check the image and fix
%it!
sizeY=fread(fileReadId,1,'*uint64');
sizeX=fread(fileReadId,1,'*uint64');
sizeT=fread(fileReadId,1,'*uint64');
if framesN==0
    framesN=sizeT;
end

switch dataType
    case 'uint8'
        dataSize=1;
    case 'uint16'
        dataSize=2;
end

sampling=fread(fileReadId,1,'*uint64');
timeStamps=zeros(framesN,1,'int64');
firstByte=30*1024+sizeX*sizeY*startT*dataSize+8*startT;
fseek(fileReadId,firstByte,-1 );

try
if strcmp(type,'frame')
    data=zeros(sizeX,sizeY,framesN,dataType);
    for t=1:1:framesN
        timeStamps(t)=fread(fileReadId,1,'*uint64');
        data(:,:,t)=fread(fileReadId,[sizeX,sizeY],['*',dataType]);
        counter=t;
    end
elseif strcmp(type,'2dROI')
    data=zeros(length(ROI(1,1):1:ROI(1,2)),length(ROI(2,1):1:ROI(2,2)),framesN,dataType);
    for t=1:1:framesN
        timeStamps(t)=fread(fileReadId,1,'*uint64');
        frame=fread(fileReadId,[sizeX,sizeY],['*',dataType]);
        data(:,:,t)=frame(ROI(1,1):1:ROI(1,2),ROI(2,1):1:ROI(2,2));
        counter=t;
    end
elseif strcmp(type,'1d')
    data=zeros(length(ROI),framesN,dataType);
    for t=1:1:framesN
        timeStamps(t)=fread(fileReadId,1,'*uint64');
        frame=fread(fileReadId,[sizeX,sizeY],['*',dataType]);
        data(:,t)=frame(ROI);
        counter=t;
    end
end
catch
    sizeT=counter;
    data(:,:,counter:end)=[];
end
fclose(fileReadId);

end