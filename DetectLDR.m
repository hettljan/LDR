function [fe,ModeShape]=DetectLDR(Matrix,fs,stopFreq,varargin)
% this function detects LDR frequency based on the pattern recognition in
% FFT slices (mode shapes), it searches for the pattern with one dominant
% vibrating region and highest max(defect)/mean(background) value
% INPUT:
%   Matrix      -   3D array with the time domain signals in 3rd dimension
%   fs          -   Sampling frequency [Hz]
%   stopFreq    -   Maximum frequency that is taken into account in the
%                   search for LDR [kHz] 
% OPTIONAL:
%   dx -            Spatial step in x direction [mm]
%   dy -            Spatial step in y direction [mm]
%   freqUnits    -   Frequency units for displaying the results
% OUTPUT:
%   fe          -   LDR frequency [kHz]
%   ModeShape   -   LDR mode shape, 2D matrix

%% INPUT PARSING
numvarargs = length(varargin);
if numvarargs > 3
    error('myfuns:somefun2Alt:TooManyInputs', ...
        'requires at most 3 optional inputs');
end
optargs = {1, 1, 'kHz'};
optargs(1:numvarargs) = varargin;
[dx, dy, freqUnits] = optargs{:};

%% PUBLIC ANNOUNCEMENT
fprintf('\n### Entering LDR detection algorithm ###\n');

%% CALCULATE THE SPATIALLY AVERAGED SPECTRA 
[Freq,FFTMatrix,~]=FFT(Matrix,fs,'lin',0,0,freqUnits);
fprintf('\nSmoothing the 2D data ...\n');
parfor i=1:size(FFTMatrix,3)
    FFTMatrix(:,:,i)=medfilt2(FFTMatrix(:,:,i),[3 3]);
end
df=Freq(2)-Freq(1);                 % frequency step in kHz
stopFreqInd=round(stopFreq/df);     % index of the stop frequency

%% LDR SEARCH ALGORITHM
processing='Im';
PeakMean=nan(stopFreqInd,1);
RestMean=nan(stopFreqInd,1);
NormFFTMatrix=nan(size(FFTMatrix,1),size(FFTMatrix,2),stopFreqInd);
for i=1:stopFreqInd
    Slice2D=FFTMatrix(:,:,i);
    switch processing
        case 'Im'
            [Stats,circle]=CircParametersFromMatrix(Slice2D,0,'auto');   % look for a reasonable circle               
            if size(Stats,1)<=2 && circle.Eccentricity <= 0.65 && circle.EquivDiameter < 25 &&  circle.EquivDiameter>4 % if it looks like circle
                Mask=zeros(size(Slice2D));    % prepare a 2D binary mask 
                x=round(circle.BoundingBox(1));
                y=round(circle.BoundingBox(2));
                w=circle.BoundingBox(3);
                h=circle.BoundingBox(4);
                if x+w <=size(Slice2D,2) && y+h <=size(Slice2D,1)
                    Mask(y:y+h,x:x+w)=1;
                elseif x+w >size(Slice2D,2) && y+h <=size(Slice2D,1)
                    Mask(y:y+h,x:end)=1;
                elseif x+w <=size(Slice2D,2) && y+h >size(Slice2D,1)
                    Mask(y:end,x:x+w)=1;
                else
                    Mask(y:end,x:end)=1;
                end
                MaskNeg=1-Mask;
                PeakMean(i)=max(Slice2D(Mask==1));
                RestMean(i)=mean(Slice2D(MaskNeg==1));       
                NormFFTMatrix(:,:,i)=Slice2D/maxND(Slice2D);
            end
        case 'Max'
            PeakMean(i)=maxND(Slice2D);
            RestMean(i)=mean(Slice2D(:));       
            NormFFTMatrix(:,:,i)=Slice2D/maxND(Slice2D);
    end
end
PNRs=10*log10(PeakMean./RestMean);      % PNR = vibration maximum divided by mean of the background 
[Gain,LDRindex]=max(PNRs);              % search for the maximum of PNR
fe=Freq(LDRindex);                      % LDR frequency
ModeShape=NormFFTMatrix(:,:,LDRindex);  % LDR mode shap
fprintf('\nLDR frequency = %.2f kHz\n',fe);
fprintf('\nPNR = %.2f dB\n',Gain);

%% ESTIMATE THE SIZE OF THE DEFECT
% % dry run to get the location of the defect
[~,circleEst]=CircParametersFromMatrix(ModeShape,0,'auto'); 
imThreshold=db2mag(-6)*maxND(ModeShape);   % define the threshold
BW = im2bw(ModeShape,imThreshold);         % convert to binary image
Stats = regionprops(BW,'EquivDiameter','Eccentricity','Extent',...
        'BoundingBox','Centroid');
for i=1:length(Stats)
    circle=Stats(i);
    if norm(circle.Centroid-circleEst.Centroid)<=5  % if the position matches the template
        fprintf('\nRadius = %.2f mm\n',circle.EquivDiameter/2)
        fprintf('\nBounding box size %.2f x %.2f mm\n',circle.BoundingBox(3),...
            circle.BoundingBox(4))
        break
    end
end

%% PLOT PNR GRAPH AND LDR MODE SHAPE
figure(1)
subplot(2,2,1)
plot(Freq(1:stopFreqInd),PNRs,'bx','MarkerSize',8);
xlabel(sprintf('Frequency [%s]',freqUnits),'FontSize',14)
ylabel('PNR [dB]','FontSize',14)
ylim([0 20])
xlim([0 stopFreq])
hold on
plot(fe,Gain,'ro','MarkerSize',12,'LineWidth',2);

% figure(2);
subplot(2,2,2)
X=(0:size(ModeShape,2)-1)*dx;
Y=(0:size(ModeShape,1)-1)*dy;
imagesc(X,Y,ModeShape);
shading interp
axis xy
axis equal
xlim([min(X) max(X)]);
ylim([min(Y) max(Y)]);
box on
hold on
rectangle('Position',circle.BoundingBox-[1 1 0 0],'LineWidth',3)
xlabel('x [mm]','FontSize',16);
ylabel('y [mm]','FontSize',16);
set(gca,'FontSize',16);
colorbar

% figure(3)
subplot(2,2,[3 4])
surf(X,Y,ModeShape,'EdgeColor','None');
shading interp
% axis equal
axis tight
% box on
xlabel('x [mm]','FontSize',16);
ylabel('y [mm]','FontSize',16);
zlabel('Norm. amplitude [-]','FontSize',16);
set(gca,'FontSize',16);
colorbar

%% INTERACTIVE PLOTTING
% PlotLDRModes(NormFFTMatrix,PNRs,1:stopFreqInd,Freq,PNRs);