close all;

savepath = '../../data/26_TU_Berlin_Mozart_Quartett_Konzerthaus/Audio/';
Nmic = 2;
Nsrc = 2;

% folder = 'Mic Distance/';
% exten = '_dist=0.3';

% folder = 'Num sources/';
% exten = '_nsrc=4';


% folder = 'Num mics/';
% exten = '_nmics=2';

folder = 'RT60/';
exten = '_rt60=5';


files = dir([savepath,folder]);
audioFiles = {};
m = 1;

for k = 1:length(files)
    if ~files(k).isdir && contains(files(k).name, [exten,'.wav'])
        audioFiles{m} = cellstr(files(k).name);
        m = m+1;
    end
end

xmic = zeros(336001,Nmic);
for k = 1:length(audioFiles)
    af = string(audioFiles{k});
    undersc = strfind(af , '_');
    nmic = str2num(extractBetween(af, undersc(1)+2, undersc(1)+2));
    nsrc = str2num(extractBetween(af, undersc(1)+1, undersc(1)+1));
    if nsrc <= Nsrc && nmic <= Nmic
        [temp,fs] = audioread(strcat(strcat(savepath, folder), af));
        xmic(:,nmic) = xmic(:,nmic) + temp;
    end
end
    
        
for nmic = 1:Nmic
    temp = xmic(:,nmic)./max(abs(xmic(:,nmic)));
    audiowrite([savepath, folder, 'mic=',num2str(nmic),exten,'.wav'], temp,fs);
end

            
    
        

        
        
        