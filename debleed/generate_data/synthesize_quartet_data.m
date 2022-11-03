%%

close all, clear all;
addpath(genpath('../../room_acoustics_simulation/RIR-Generator/.'));   %rir toolbox
soundpath = '../../data/26_TU_Berlin_Mozart_Quartett_Konzerthaus/Audio/';  %anechoic string quartet recordings
savemat = '../../data/26_TU_Berlin_Mozart_Quartett_Konzerthaus/Mat files/';   %save generated mat files
fname = {'Va','Vcl','Vl1','Vl2'};
fs = 48000;
ch = 'mic';     %which scenario do I want to simulate?


%% RIR config, everything is in SI units

Nmic = 4;  %number of mics
c = 340;    %sound velocity
x = 3;  %length
y = 4;  %width
z = 3.25;  %height
L = [x y z];    %room dimensions (x,y,z)

% place quartet in a semicircle
xcen = 1.5;
ycen = 1; %2.5;
rad = 0.5; %0.4;
theta = linspace(0,pi,100);
zcen = 1;

angpos = (linspace(0, pi ,Nmic) + 0.05).';
xpos = xcen+ rad*cos(angpos);
ypos = ycen + rad*sin(angpos);
zpos = zcen*ones(Nmic,1);
zpos(2) = zcen - 0.2; %cello (should be lower than others)
s = [xpos ypos zpos];   %source positions (x,y,z)
r = s - repmat([0 0.2 0], Nmic, 1);   %receiver positions for close mics (x,y,z)

beta = 0.8;     %RT60(s)
n = 2^8;   %number of samples of RIR
mtype = 'omnidirectional';  %mic type
order = 3;  %reflection order
dim = 3;    %room dimensions
orientation = [0 0];    %direction where mic is pointed (azimuth and elevation angle)
hp_filter = 1;  %enable high pass filter


% use a different portion of the sound file for calibration
calib_start = round(4*fs);
calib_stop = round(6.5*fs);


%% choose which scenario to simulate

switch (ch)
    
   %% vary source-mic distance
    case 'dist'
        
        dis = 0.1:0.1:0.5;    %receiver array
        N = 2;
        start = round(35*fs);
        stop = round(42*fs);
        sim.dis = 0;
        sim.N = N;
        sim.src_name = ['Va','Vcl','Vl1,','Vl2'];
        sim.h_ideal = zeros(n,N,N);
        sim.src = struct([]);
        sim.calib = struct([]);
        sim.mic = zeros(stop-start+1,N);
        Sim = repmat(sim,1,length(dis));
        savepath = [soundpath,'Mic Distance/'];

        for k = 1:length(dis)

            disp(['Generating data for src-mic distance ',num2str(dis(k)),' m']);
            r = s - repmat([0 dis(k) 0], Nmic, 1); 
            Sim(k).dis = dis(k);

            %% plot mic and source configureation
            figure(1);
            sp = scatter3(s(:,1), s(:,2), s(:,3));hold on;grid on;
            sp.MarkerEdgeColor = 'b';
            rp = scatter3(r(:,1), r(:,2), r(:,3), 'x');hold on;
            rp.MarkerEdgeColor = 'r';
            plot3(xcen+rad*cos(theta), ycen+rad*sin(theta),zcen*ones(length(theta),1),'-k','LineWidth',1.2);hold on;
            xlim([0 L(1)]);ylim([0 L(2)]);zlim([0 L(3)]);
            drawnow;
            title('Room configuration');
            legend('Sources','Mics');


            %% loop to create RIRs and convolve them with audio

            %loop over sources
            for i = 1:N
                [x,fs] = audioread([soundpath,fname{i},'.wav']);
                x = x(:,1); %convert stereo to mono
                x_trim = trim_audio(x,fs,start,stop);
                x_calib = trim_audio(x,fs,calib_start, calib_stop);


                h = rir_generator(c, fs, r(1:N,:), s(i,:), L, beta, n, mtype, order, dim, orientation, hp_filter);
                figure(2); plot(h.');title(['Distance =' ,num2str(dis(k))]);
                Sim(k).src{i}.clean = x_trim;
                Sim(k).src{i}.mic = zeros(length(x_trim),N);
                %loop over receivers
                for j = 1:N
                     y = fftfilt(h(j,:),x_trim);
                     Sim(k).src{i}.mic(:,j) = y;
                     Sim(k).calib(1).src{i}.mic(:,j) = fftfilt(h(j,:), x_calib);
                     Sim(k).h_ideal(:,j,i) = h(j,:);
                     audiowrite([savepath,fname{i},'_',num2str(i),num2str(j),'_mic_dist=',num2str(dis(k)),'.wav'],y,fs);
                end 
            end

            for i = 1:N
                for j = 1:N
                m_trim = Sim(k).src{j}.mic(:,i);
                Sim(k).mic(:,i) = Sim(k).mic(:,i) + m_trim;
                end
            end
        end

        save([savemat, 'quartet_synthesized_data_var_distance_Nsrc=',num2str(N),'.mat'],'Sim');



    %% vary volume and RT60
    case 'vol'
        dis = 0.2;
        r = s - repmat([0 dis 0], Nmic, 1);
        V = 2*logspace(1,2,5);  %logarithmically placed volume
        x = nthroot(V,3).'; y = x; z = x;
        beta = [0.5 1 2 5 8];   %increasing RT60
        n = 2^10;    %use at least 5ms of data for convolution

        N = 2;
        start = round(35*fs);
        stop = round(42*fs);
        sim.V = 0;
        sim.beta = 0;
        sim.N = N;
        sim.src_name = ['Va','Vcl','Vl1,','Vl2'];
        sim.h_ideal = zeros(n,N,N);
        sim.calib = struct([]);
        sim.src = struct([]);
        sim.mic = zeros(stop-start+1,N);
        Sim = repmat(sim,1,length(beta));
        savepath = [soundpath,'RT60/'];
        

        for k = 1:length(V)

            disp(['Generating data Volume = ',num2str(V(k)),' m^3, RT60 = ', num2str(beta(k)), ' s']);
            L = round([x(k), y(k), z(k)]*100)/100;
            Sim(k).V = prod(L);
            Sim(k).beta = beta(k);
            
            %% plot mic and source configureation
            figure(1);
            sp = scatter3(s(:,1), s(:,2), s(:,3));hold on;grid on;
            sp.MarkerEdgeColor = 'b';
            rp = scatter3(r(:,1), r(:,2), r(:,3), 'x');hold on;
            rp.MarkerEdgeColor = 'r';
            plot3(xcen+rad*cos(theta), ycen+rad*sin(theta),zcen*ones(length(theta),1),'-k','LineWidth',1.2);hold on;
            xlim([0 L(1)]);ylim([0 L(2)]);zlim([0 L(3)]);
            drawnow;
            title('Room configuration');
            legend('Sources','Mics');

            %% loop to create RIRs and convolve them with audio
            for i = 1:N
                [xmic,fs] = audioread([soundpath,fname{i},'.wav']);
                xmic = xmic(:,1); %convert stereo to mono
                x_trim = trim_audio(xmic,fs,start,stop);
                x_calib = trim_audio(xmic,fs,calib_start, calib_stop);


                h = rir_generator(c, fs, r(1:N,:), s(i,:), L, beta(k), n, mtype, order, dim, orientation, hp_filter);
                Sim(k).src{i}.clean = x_trim;
                Sim(k).src{i}.mic = zeros(length(x_trim),N);
                
                %loop over receivers
                for j = 1:N
                     yfilt = fftfilt(h(j,:),x_trim);
                     Sim(k).src{i}.mic(:,j) = yfilt;
                     Sim(k).h_ideal(:,j,i) = h(j,:);
                     Sim(k).calib(1).src{i}.mic(:,j) = fftfilt(h(j,:), x_calib);
                     audiowrite([savepath,fname{i},'_',num2str(i),num2str(j),'_rt60=',num2str(beta(k)),'.wav'],yfilt,fs);
                end 
            end

            for i = 1:N
                for j = 1:N
                m_trim = Sim(k).src{j}.mic(:,i);
                Sim(k).mic(:,i) = Sim(k).mic(:,i) + m_trim;
                end
            end
        end

        save([savemat, 'quartet_synthesized_data_var_rt60+vol_Nsrc=',num2str(N),'.mat'],'Sim');


    %% vary rt60
    case 'rt60'
        dis = 0.2;
        r = s - repmat([0 dis 0], Nmic, 1); 
        beta = [0.5 1 2 5 8];
        N = 2;
        start = round(35*fs);
        stop = round(42*fs);
        sim.rt60 = 0;
        sim.N = N;
        sim.src_name = ['Va','Vcl','Vl1,','Vl2'];
        sim.h_ideal = zeros(n,N,N);
        sim.src = struct([]);
        sim.mic = zeros(stop-start+1,N);
        Sim = repmat(sim,1,length(beta));
        savepath = [soundpath,'RT60/'];

        for k = 1:length(beta)

            disp(['Generating data RT60 = ',num2str(beta(k)),' s']);
            Sim(k).rt60 = beta(k);

            %% loop to create RIRs and convolve them with audio
            for i = 1:N
                [x,fs] = audioread([soundpath,fname{i},'.wav']);
                x = x(:,1); %convert stereo to mono
                x_trim = trim_audio(x,fs,start,stop);

                h = rir_generator(c, fs, r(1:N,:), s(i,:), L, beta(k), n, mtype, order, dim, orientation, hp_filter);
                Sim(k).src{i}.clean = x_trim;
                Sim(k).src{i}.mic = zeros(length(x_trim),N);
                %loop over receivers
                for j = 1:N
                     y = fftfilt(h(j,:),x_trim);
                     Sim(k).src{i}.mic(:,j) = y;
                     Sim(k).h_ideal(:,j,i) = h(j,:);
                     audiowrite([savepath,fname{i},'_',num2str(i),num2str(j),'_rt60=',num2str(beta(k)),'.wav'],y,fs);
                end 
            end

            for i = 1:N
                for j = 1:N
                m_trim = Sim(k).src{j}.mic(:,i);
                Sim(k).mic(:,i) = Sim(k).mic(:,i) + m_trim;
                end
            end
        end

        save([savemat, 'quartet_synthesized_data_var_rt60_Nsrc=',num2str(N),'.mat'],'Sim');



    %% vary number of microphones
    case 'mic'
        dis = 0.2;
        Nsrc = 2;
        s = s(1:Nsrc,:);
        Nmic = 1:4;
        start = round(35*fs);
        stop = round(42*fs);

        sim.Nsrc = Nsrc;
        sim.Nmic = zeros(length(Nmic),1);
        sim.src_name = ['Va','Vcl','Vl1,','Vl2'];

        sim.src = struct([]);
        sim.calib = struct([]);
        sim.mic = [];
        sim.h_ideal = [];
        Sim = repmat(sim,1,length(Nmic));
        savepath = [soundpath,'Num Mics/'];
        ang = linspace(0,pi/2,max(Nmic)).';

        for k = 1:length(Nmic)

            disp(['Generating data for number of mics ',num2str(Nmic(k))]);

            Sim(k).Nmic = Nmic(k);
            mang = ang(1:k);   % angle of mic in xy plane
            ss = repelem(s,Nmic(k),1);
            rr = repmat([dis*sin(mang) dis*cos(mang) zeros(length(mang),1)], Nsrc,1);
            r = ss - rr;
            cart_dist = sqrt(sum((ss - r).^2,2));

            figure(1);clf;
            sp = scatter3(s(:,1), s(:,2), s(:,3));hold on;grid on;
            sp.MarkerEdgeColor = 'b';
            rp = scatter3(r(:,1), r(:,2), r(:,3), 'x');hold on;
            rp.MarkerEdgeColor = 'r';
            plot3(xcen+rad*cos(theta), ycen+rad*sin(theta),zcen*ones(length(theta),1),'-k','LineWidth',1.2);hold off;
            xlim([0 L(1)]);ylim([0 L(2)]);zlim([0 L(3)]);
            drawnow;
            title('Room configuration');
            legend('Sources','Mics');

            Sim(k).h_ideal = zeros(n,Nmic(k)*Nsrc,Nsrc);
            Sim(k).mic = zeros(stop-start+1,Nmic(k)*Nsrc);

            %% loop to create RIRs and convolve them with audio

            for i = 1:Nsrc
                [x,fs] = audioread([soundpath,fname{i},'.wav']);
                x = x(:,1); %convert stereo to mono
                x_trim = trim_audio(x,fs,start,stop);
                x_calib = trim_audio(x,fs,calib_start, calib_stop);

                h = rir_generator(c, fs, r, s(i,:), L, beta, n, mtype, order, dim, orientation, hp_filter);
                Sim(k).src{i}.clean = x_trim;
                Sim(k).src{i}.mic = zeros(length(x_trim),Nmic(k));
                %loop over receivers
                for j = 1:Nmic(k)*Nsrc
                     y = fftfilt(h(j,:),x_trim);
                     Sim(k).src{i}.mic(:,j) = y;
                     Sim(k).calib(1).src{i}.mic(:,j) = fftfilt(h(j,:), x_calib);
                     Sim(k).h_ideal(:,j,i) = h(j,:);
                     audiowrite([savepath,fname{i},'_',num2str(i),num2str(j),'_nmics=',num2str(Nmic(k)),'.wav'],y,fs);
                end 
            end

            for i = 1:Nmic(k)*Nsrc
                for j = 1:Nsrc
                    m_trim = Sim(k).src{j}.mic(:,i);
                    Sim(k).mic(:,i) = Sim(k).mic(:,i) + m_trim;           
                end

            end
        end

        save([savemat, 'quartet_synthesized_data_var_mics_Nsrc=', num2str(Nsrc),'.mat'],'Sim');


    %% vary number of sources
    case 'src'

        dis = 0.2;
        Nsrc = 2:4;
        Nmic = Nsrc;    %determined system
        start = round(35*fs);
        stop = round(42*fs);

        sim.Nsrc = zeros(length(Nsrc),1);
        sim.Nmic = zeros(length(Nsrc),1);
        sim.src_name = ['Va','Vcl','Vl1,','Vl2'];

        sim.src = struct([]);
        sim.mic = [];
        sim.h_ideal = [];
        sim.calib = struct([]);
        Sim = repmat(sim,1,length(Nsrc));
        savepath = [soundpath,'Num sources/'];

        for k = 1:length(Nsrc)

            disp(['Generating data for number of sources ',num2str(Nsrc(k))]);

            Sim(k).Nsrc = Nsrc(k);
            Sim(k).Nmic = Nmic(k);
            ss = s(1:Nsrc(k),:);
            rr = ss - repmat([0 dis 0], Nmic(k), 1);   %receiver positions for close mics (x,y,z)

            figure(1);clf;
            sp = scatter3(ss(:,1), ss(:,2), ss(:,3));hold on;grid on;
            sp.MarkerEdgeColor = 'b';
            rp = scatter3(rr(:,1), rr(:,2), rr(:,3), 'x');hold on;
            rp.MarkerEdgeColor = 'r';
            plot3(xcen+rad*cos(theta), ycen+rad*sin(theta),zcen*ones(length(theta),1),'-k','LineWidth',1.2);hold off;
            xlim([0 L(1)]);ylim([0 L(2)]);zlim([0 L(3)]);
            drawnow;
            title('Room configuration');
            legend('Sources','Mics');

            Sim(k).h_ideal = zeros(n,Nmic(k),Nsrc(k));
            Sim(k).mic = zeros(stop-start+1,Nmic(k));

            %% loop to create RIRs and convolve them with audio

            for i = 1:Nsrc(k)
                [x,fs] = audioread([soundpath,fname{i},'.wav']);
                x = x(:,1); %convert stereo to mono
                x_trim = trim_audio(x,fs,start,stop);
                x_calib = trim_audio(x,fs,calib_start, calib_stop);

                h = rir_generator(c, fs, rr, ss(i,:), L, beta, n, mtype, order, dim, orientation, hp_filter);
                Sim(k).src{i}.clean = x_trim;
                Sim(k).src{i}.mic = zeros(length(x_trim),Nmic(k));

                %loop over receivers
                for j = 1:Nmic(k)
                     y = fftfilt(h(j,:),x_trim);
                     Sim(k).src{i}.mic(:,j) = y;
                     Sim(k).calib(1).src{i}.mic(:,j) = fftfilt(h(j,:), x_calib);
                     Sim(k).h_ideal(:,j,i) = h(j,:);
                     audiowrite([savepath,fname{i},'_',num2str(i),num2str(j),'_nsrc=',num2str(Nsrc(k)),'.wav'],y,fs);
                end 
            end

            for i = 1:Nmic(k)
                for j = 1:Nsrc(k)
                    m_trim = Sim(k).src{j}.mic(:,i);
                    Sim(k).mic(:,i) = Sim(k).mic(:,i) + m_trim;           
                end
            end
        end

        save([savemat, 'quartet_synthesized_data_var_src.mat'],'Sim');
    
    %% if nothing matches
    otherwise
        error('Wrong choice!');
        
end