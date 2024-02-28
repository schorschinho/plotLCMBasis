% plotLCMBasis.m
% Georg Oeltzschner, Johns Hopkins University 2024.
%
% USAGE:
% out = plotLCMBasis(basisSetFile, stagFlag, ppmmin, ppmmax, xlab, ylab, figTitle)
% 
% DESCRIPTION:
% Creates a figure showing all a basis functions in LCModel basis set.
% 
% INPUTS:
% basisSetFile  = LCModel basis set (.basis)
% stagFlag  = flag to decide whether basis functions should be plotted
%               vertically staggered or simply over one another (optional.
%               Default: 1 = staggered; 0 = not staggered)
% ppmmin    = lower limit of ppm scale to plot (optional.  Default = 0.2 ppm).
% ppmmax    = upper limit of ppm scale to plot (optional.  Default = 5.2 ppm).
% xlab      = Label for the x-axis (optional.  Default = 'Chemical shift (ppm)');
% ylab      = label for the y-axis (optional.  Default = '');
% figTitle  = label for the title of the plot (optional.  Default = '');

function out = plotLCMBasis(basisSetFile, stagFlag, ppmmin, ppmmax, xlab, ylab, figTitle)

% Parse input arguments
if nargin<7
    figTitle = '';
    if nargin<6
        ylab='';
        if nargin<5
            xlab='Chemical shift (ppm)';
            if nargin<4
                ppmmax=5.2;
                if nargin<3
                    ppmmin=0.2;
                    if nargin<2
                        stagFlag = 1;
                        if nargin<1
                            error('ERROR: no input basis set specified.  Aborting!!');
                        end
                    end
                end
            end
        end
    end
end

% Generate a new figure and keep the handle memorized
out = figure;

% Read in the basis set file
basisSet = io_readlcmraw_basis(basisSetFile);

% Reorder alphabetically (ignoring case)
[~, newOrder]   = sort(lower(fieldnames(basisSet)));
basisSet        = orderfields(basisSet, newOrder);

% Define number of basis functions
metNamesArray   = fieldnames(basisSet);
nBasisFct       = length(metNamesArray);
nPts            = basisSet.(metNamesArray{1}).sz(1);

% Extract FIDs & metadata
fidsArray       = zeros(nPts, nBasisFct);
for kk = 1:nBasisFct
    fidsArray(:,kk) = basisSet.(metNamesArray{kk}).fids;
end
specsArray      = fftshift(fft(fidsArray,[],1),1);
dwellTime       = basisSet.(metNamesArray{1}).dwelltime;
txfrq           = basisSet.(metNamesArray{1}).txfrq;
spectralWidth   = 1/dwellTime;
f               = [(-spectralWidth/2)+(spectralWidth/(2*nPts)):spectralWidth/(nPts):(spectralWidth/2)-(spectralWidth/(2*nPts))];
ppm             = f/(txfrq*1e-6) + 4.68;

% Determine stag and plot
if stagFlag==1

    % Staggered plots will be separated by the mean of the
    % maximum across all spectra
    stag = mean(max(real(specsArray)));
    
    % Loop over all basis functions
    hold on
    for kk = 1:nBasisFct
        plot(ppm, real(specsArray(:,kk) - kk*stag), 'k');
        % Instead of a MATLAB legend, annotate each line separately with the
        % name of the metabolite
        text(ppmmin, - kk*stag, metNamesArray{kk}, 'FontSize', 14);
    end
    hold off
    
    set(gca, 'XDir', 'reverse', 'XLim', [ppmmin, ppmmax], 'YLim', [-stag*(nBasisFct+1), max(real(specsArray(:,1)))]);

else

    % Loop over all basis functions
    hold on
    for kk = 1:nBasisFct
        plot(ppm, real(specsArray(:,kk)));
    end
    legend(metNamesArray);
    hold off

    set(gca, 'XDir', 'reverse', 'XLim', [ppmmin, ppmmax]);

end

% Common style for all outputs
set(gca, 'LineWidth', 1, 'TickDir', 'out');
set(gca, 'FontSize', 16);
% If no y caption, remove y axis
if isempty(ylab)
    set(gca, 'YColor', 'w');
else
    set(gca, 'YColor', 'k');
end
% Black axes, white background
set(gca, 'XColor', 'k');
set(gca, 'Color', 'w');
set(gcf, 'Color', 'w');
box off;
title(figTitle);
xlabel(xlab, 'FontSize', 20);
ylabel(ylab, 'FontSize', 20);

% Set linewidth coherently
Fig1Ax1 = get(out, 'Children');
Fig1Ax1Line1 = get(Fig1Ax1, 'Children');
if iscell(Fig1Ax1Line1)
    Fig1Ax1Line1 = Fig1Ax1Line1(~cellfun('isempty', Fig1Ax1Line1));
    Fig1Ax1Line1 = Fig1Ax1Line1{1};
end
set(Fig1Ax1Line1, 'LineWidth', 1);


end



% io_readlcmraw_basis.m
% Jamie Near, McGill University 2014.
%
% USAGE:
% out=io_readlcmraw_basis(filename);
%
% DESCRIPTION:
% Reads entire LCModel .basis file into multiple FID-A data structures in MATLAB.
%
% INPUTS:
% filename   = filename of LCModel .basis file.
% conjugate       = apply complex conjugate
%
% OUTPUTS:
% out        = Input basis set saved as a structure in which each field is
%               an individual metabolite basis spectrum in FID-A structure 
%               format.

function [out]=io_readlcmraw_basis(filename,conjugate)

if nargin < 2
    conjugate = 1;
end

% Begin to read the data.
fid=fopen(filename);
line=fgets(fid);

%look for FWHMBA
fwhmba_index=contains(line,'FWHMBA');
while ~fwhmba_index;
    line=fgets(fid);
    fwhmba_index=contains(line,'FWHMBA');
    line = GetNumFromString(line);
end
linewidth=str2num(line);

%look for HZPPM
hzpppm_index=contains(line,'HZPPPM');
while ~hzpppm_index;
    line=fgets(fid);
    hzpppm_index=contains(line,'HZPPPM');
    line = GetNumFromString(line);
end
hzpppm=str2num(line);
Bo=hzpppm/42.577;
linewidth=linewidth*hzpppm;

%look for TE
te_index=contains(line,'ECHOT');
while ~te_index;
    line=fgets(fid);
    te_index=contains(line,'ECHOT');
    line = GetNumFromString(line);
end
te=str2num(line);

%look for spectral width
badelt_index=contains(line,'BADELT');
while ~badelt_index;
    line=fgets(fid);
    badelt_index=contains(line,'BADELT');
    line = GetNumFromString(line);
end
dwelltime=str2num(line);
spectralwidth=1/dwelltime;
fileEnd=false;

while ~feof(fid)
     %look for a center frequency
    ppmsep_index=contains(line,'PPMSEP');    
    while ~ppmsep_index;
        line=fgets(fid);
        if ischar(line)
            ppmsep_index=contains(line,'PPMSEP');

            % ARC 2023-06-14 'METABO' also matches things like "METABO_CONTAM"
            % and "METABO_SINGLET", which breaks eval'd code later. Hence,
            % catch 'METABO' then reject 'METABO_':
            if contains(line,'METABO') && ~contains(line, 'METABO_')
                break
            end
        end       
    end
    if ppmsep_index
        line = GetNumFromString(line);
        centerFreq=str2num(line);
    else
        centerFreq = [];
    end
    
    %Look for the metabolite name
    metab_index=contains(line,'METABO') && ~contains(line, 'METABO_'); % ARC 2023-06-14
    while ~metab_index;
        line=fgets(fid);
        if ischar(line)
            metab_index=contains(line,'METABO') && ~contains(line, 'METABO_'); % ARC 2023-06-14
            if metab_index
                break
            end
        end
    end

    % ARC 2023-06-14 : regular expression to extract name part (may need to be adapted in case of odd characters in the metabolite name)
    pat=regexp(line,'^\s*METABO\s*=\s*[''"]([-_+A-Za-z0-9]+)\s*[''"].*$','tokens');

    if isempty(pat)
        continue;  % failed to extract a meaningful name... move on to the next metab
    end

    metabName=pat{1}{1};
    
    hdrEnd_index=contains(line,'$END');
    while ~hdrEnd_index;
        line=fgets(fid);
        if ischar(line)
            hdrEnd_index=contains(line,'$END');
            if hdrEnd_index
                break
            end
        end
    end
    
    line=fgets(fid);
    
    nmused_index=contains(line,'$NMUSED');
    basis_index=contains(line,'$BASIS');
    linenum=1;
    RF=[];
    % If the line is empty skip it
    while ~isempty(line) && ~nmused_index && ~basis_index && ~fileEnd
        %dataline=line(1:semicol_index-2);
        [A,count, errmsg, nextindex] = sscanf(line, '%f', inf);
        % If read failed, output the error
        if ~isempty(errmsg);
            fclose(fid);
            error('READLCMRAW_BASIS failed with read error: %s', errmsg);
        end
        % Store the read values into rf array
        RF = [ RF ; A ];
        if feof(fid)
            fileEnd=true;
        end
        linenum = linenum + 1;
        line=fgets(fid);
        if ischar(line)
            nmused_index=contains(line,'$NMUSED');
            basis_index=contains(line,'$BASIS');
        end
    end
    specs=RF(1:2:end) + 1i*RF(2:2:end);

    % GO 2022/01/24
    % LCModel uses negative BADELT values to encrypt basis sets
    % (LCModel.f source code lines 3666 and following)
    if dwelltime < 0
        dix = 1499;
        for rr = 1:length(specs)
            [randomresult, dix] = lcmodelrng(dix);
            specs(rr) = -specs(rr) .* exp(-20*randomresult + 10);
        end
    end

    if conjugate
        specs=flipud(fftshift(conj(specs),1));
    else
        specs=(fftshift(specs,1));
    end
    vectorsize=length(specs);
    sz=[vectorsize 1];
    if mod(vectorsize,2)==0
        %disp('Length of vector is even.  Doing normal conversion');
        fids=ifft(ifftshift(specs,1),[],1);
    else
        %disp('Length of vector is odd.  Doing circshift by 1');
        fids=ifft(circshift(ifftshift(specs,1),1),[],1);
    end
    f=[(-spectralwidth/2)+(spectralwidth/(2*sz(1))):spectralwidth/(sz(1)):(spectralwidth/2)-(spectralwidth/(2*sz(1)))];
    ppm=f/(Bo*42.577);
    ppm=ppm+4.68;
    t=[dwelltime:dwelltime:vectorsize*dwelltime];
    txfrq=hzpppm*1e6;
    metabName = RemoveWhiteSpaces(metabName);
    if strcmp(metabName,'-CrCH2')
        metabName = 'CrCH2';
    end
    if strcmp(metabName,'2HG')
        metabName = 'bHG';
    end
    eval(['out.' metabName '.fids=fids;']);
    eval(['out.' metabName '.specs=specs;']);
    eval(['out.' metabName '.sz=[vectorsize 1 1 1];']);
    eval(['out.' metabName '.n=vectorsize;']);
    eval(['out.' metabName '.spectralwidth=abs(spectralwidth);']);
    eval(['out.' metabName '.sz=sz;']);
    eval(['out.' metabName '.Bo=Bo;']);
    eval(['out.' metabName '.te=te;']);
    eval(['out.' metabName '.tr=[];']);
    eval(['out.' metabName '.dwelltime=abs(1/spectralwidth);']);
    eval(['out.' metabName '.linewidth=linewidth;']);
    eval(['out.' metabName '.ppm=ppm;']);
    eval(['out.' metabName '.t=t;']);
    eval(['out.' metabName '.txfrq=txfrq;']);
    if ~isempty(centerFreq)
            eval(['out.' metabName '.centerFreq=centerFreq;']);
    end
    eval(['out.' metabName '.date=date;']);
    eval(['out.' metabName '.seq='''';']);
    eval(['out.' metabName '.sim='''';']);
    eval(['out.' metabName '.dims.t=1;']);
    eval(['out.' metabName '.dims.coils=0;']);
    eval(['out.' metabName '.dims.averages=0;']);
    eval(['out.' metabName '.dims.subSpecs=0;']);
    eval(['out.' metabName '.dims.extras=0;']);
    eval(['out.' metabName '.averages=1;']);
    eval(['out.' metabName '.flags.writtentostruct=1;']);
    eval(['out.' metabName '.flags.gotparams=1;']);
    eval(['out.' metabName '.flags.leftshifted=1;']);
    eval(['out.' metabName '.flags.filtered=0;']);
    eval(['out.' metabName '.flags.zeropadded=0;']);
    eval(['out.' metabName '.flags.freqcorrected=0;']);
    eval(['out.' metabName '.flags.phasecorrected=0;']);
    eval(['out.' metabName '.flags.averaged=1;']);
    eval(['out.' metabName '.flags.addedrcvrs=1;']);
    eval(['out.' metabName '.flags.subtracted=1;']);
    eval(['out.' metabName '.flags.writtentotext=1;']);
    eval(['out.' metabName '.flags.downsampled=0;']);
    eval(['out.' metabName '.flags.isFourSteps=0;']);
    
end

fclose(fid);
end

% RF=RF';
% rf(:,1)=RF(:,2)*180/pi;
% rf(:,2)=RF(:,1);
% rf(:,3)=ones(length(RF(:,1)),1);

function str3 = GetNumFromString(str)
str1 = regexprep(str,'[,;=]', ' ');
str2 = regexprep(regexprep(str1,'[^- 0-9.eE(,)/]',''), ' \D* ',' ');
str3 = regexprep(str2, {'\.\s','\E\s','\e\s','\s\E','\s\e'},' ');
end

function str = RemoveWhiteSpaces(str)
    pattern = '[ \t\n]'; % Match zero or more spaces, tabs, or newlines, followed by a double quote
    replacement = ''; % Replace the matched string with just a double quote
    str = regexprep(str, pattern, replacement);
end


% GO 2022/01/24
% LCModel uses negative BADELT values to encrypt basis sets
% (LCModel.f source code lines 3666 and following)
%++++++++++++++++ doubleprecision VERSION 2DP (MAR 1984) ++++++++++++++    4222
%  FUNCTION RANDOM.  PRODUCES A PSEUDORANDOM REAL ON THE OPEN INTERVAL      4223
%      (0.,1.).                                                             4224
%  DIX (IN doubleprecision) MUST BE INITIALIZED TO A WHOLE NUMBER          4225
%      BETWEEN 1.0D0 AND 2147483646.0D0 BEFORE THE FIRST CALL TO RANDOM       4226
%      AND NOT CHANGED BETWEEN SUCCESSIVE CALLS TO RANDOM.                  4227
%  BASED ON L. SCHRAGE, ACM TRANS. ON MATH. SOFTWARE 5, 132 (1979).         4228
%-----------------------------------------------------------------------    4229
function [randomresult,dix]=lcmodelrng(dix);

a   = 16807.0d0;
b15 = 32768.d0;
b16 = 65536.0d0;
p   = 2147483647.0d0;

 %                                                                           4231
 %  PORTABLE RANDOM NUMBER GENERATOR                                         4232
 %   USING THE RECURSION                                                     4233
 %    DIX = DIX*A MOD P                                                      4234
 %                                                                           4235
 %                                                                           4237
 %  7**5, 2**15, 2**16, 2**31-1                                              4238
                                                             
 %  GET 15 HI ORDER BITS OF DIX                                              4241
 xhi = dix./b16;
 xhi = xhi - rem(xhi,1.0d0);
 %  GET 16 LO BITS IF DIX AND FORM LO PRODUCT                                4244
 xalo =(dix-xhi.*b16).*a;
 %  GET 15 HI ORDER BITS OF LO PRODUCT                                       4246
 leftlo = xalo./b16;
 leftlo = leftlo - rem(leftlo,1.0d0);
 %  FORM THE 31 HIGHEST BITS OF FULL PRODUCT                                 4249
 fhi = xhi.*a + leftlo;
 %  GET OVERFLO PAST 31ST BIT OF FULL PRODUCT                                4251
 k = fix(fhi./b15);
 k = fix(k - rem(k,1.0d0));
 %  ASSEMBLE ALL THE PARTS AND PRESUBTRACT P                                 4254
 %   THE PARENTHESES ARE ESSENTIAL                                           4255
 dix =(((xalo-leftlo.*b16)-p)+(fhi-k.*b15).*b16) + k;
 %  ADD P BACK IN IF NECESSARY                                               4257
 if(dix < 0.0d0)
    dix = dix + p;
 end
 %  MULTIPLY BY 1/(2**31-1)                                                  4259
 randomresult = dix.*4.656612875d-10;
end %function random