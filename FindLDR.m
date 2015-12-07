clear all

%% SCAN LOAD SETTINGS
folder='C:\Users\u0088749\KULeuven\PROJECTS\ALAMSA\Signals\FBH\LargeSample';
file='Setup6Sensor4.tdms';
fileName=fullfile(folder,file);
% chNamePattern='x %f y %f/Amp1.000V';
chNamePattern='x %f y %f/Phase0'; 
fs=2e6;              % Sampling frequency [Hz]
SpatialSteps=[1 1];   % Spatial step in [mm] 

%% SIGNAL PROCESSING SETTINGS
detrendOn=1;            % remove constant trend
filterOn=0;             % enable the filter
filterType='Band';
Fstop1=80e3;           % Low stop freqeuncy
Fpass1=100e3;           % Low pass frequency
Fpass2=200e3;           % High pass frequency
Fstop2=220e3;           % High stop freqeuncy

%% VISUALIZATION SETTINGS
freqUnit='kHz';         
medianFiltOn=1;         % allow the median filter for smoothing frames for video?
TimeLims=[nan nan];     % indeces used to crop the matrix in time
stopFreq=100;           % Upper frequency limit in the given units

%% LOAD DATA
Matrix=LoadCscanDataTDMS(fileName,SpatialSteps(1),SpatialSteps(2),chNamePattern);

%% SQUEEZE AND ROTATE
fprintf('\nSqueezing and rotating data ...\n');
Matrix=rot90(Matrix);   % rotate by 90 deg to get the correct orientation
Matrix=flipud(Matrix(1:end,:,:));

%% CROP TIME DOMAIN
if isnan(TimeLims(1)) == 0      % if the value is non-nan
    fprintf('\nCropping signals ...\n');
    if isnan(TimeLims(2)) == 0  % if the value is non-nan
        Matrix=Matrix(:,:,TimeLims(1):TimeLims(2));
    else
        Matrix=Matrix(:,:,TimeLims(1):end);
    end
end 

%% DESIGN THE FILTER
if filterOn == 1
    fprintf('\nDesigning Filter ...\n');
    switch filterType
        case 'Band'
            Hd = BandPassFIR(Fstop1,Fpass1,Fpass2,Fstop2,fs);
        case 'Low'
            Hd = LowPassFIR(Fpass2,Fstop2,fs);
    end
fvtool(Hd)    
end

%% DETREND, FILTER, NORMALIZE
[m,n,p]=size(Matrix);
Matrix=reshape(Matrix,[],p);        % flaten to 2D matrix
Matrix=Matrix';                     % transpose to columnwise vector
if detrendOn == 1
    fprintf('\nDetrend ... \n');
    Matrix=detrend(Matrix,'constant');
end
if filterOn == 1
    fprintf('\nFiltering ... \n');
    Matrix= filtfilt(Hd,Matrix);
end
Matrix=Matrix';                     % retranspose 
Matrix=reshape(Matrix,m,n,p);       % shape into th

%% IDENTIFY THE LDR MODE AND FREQUENCY
[fe,ModeShape]=DetectLDR(Matrix,fs,stopFreq,1,1);

%% 
saveas(gcf,'C:\Users\u0088749\Desktop\untitled.png');
