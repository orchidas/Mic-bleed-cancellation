close all, clear all;

savepath = '../../data/26_TU_Berlin_Mozart_Quartett_Konzerthaus/Audio/';
N = 2;
% folder = ['Num mics/Num sources/N=',num2str(N),'/'];
% folder = ['Mic Distance/Num Sources/N=',num2str(N),'/'];
% folder = ['RT60/Num sources/N=',num2str(N),'/'];
% folder = ['Num sources/'];
calib_type = 'spec-ratio';


dirs = {'Ideal','MCWF', 'GEVD-MWF', 'MLE', 'MAP'};
% dirs = {'MAP'};
hyp = [0 1 100];

for k = 1:length(dirs)
   
    if strcmp(dirs{k}, 'MLE') || strcmp(dirs{k}, 'MAP')
        str_add = ['_',calib_type,'_', lower(dirs{k})];
        for j = 1:length(hyp)
            hyp_folder = ['sigma=',num2str(hyp(j)),'/'];
            fullpath = [savepath, folder, dirs{k}, '/',hyp_folder];
            files = dir(fullpath);
            for i = 1:length(files)
                file = files(i).name;
%                 rename = [file(1:end-4), str_add, '.wav'];
                if endsWith(file, [str_add, '.wav'])
%                     movefile(fullfile(fullpath,file), fullfile(fullpath,rename));
                    [x,fs] = audioread([fullpath,file]);
                    y = x./max(abs(x));
                    audiowrite([savepath, folder, dirs{k}, '/',hyp_folder,file], y, fs);
                end
            end
        end
    else
         files = dir([savepath, folder, dirs{k}, '/']);
         for i = 1:length(files)
            file = files(i).name;
            if endsWith(file, '.wav')
                [x,fs] = audioread([savepath, folder, dirs{k}, '/',file]);
                y = x./max(abs(x));
                audiowrite([savepath, folder, dirs{k}, '/',file], y, fs);
            end
         end
    end

end
