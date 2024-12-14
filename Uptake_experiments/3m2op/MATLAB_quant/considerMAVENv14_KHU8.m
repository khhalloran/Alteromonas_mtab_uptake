function [sampleNames, keepGoodData, exportCurves] = considerMAVENv14(CSVfileFromMAVEN,sampleInfoFile,setQuality,pruneData,ionMode,SRMfile)
%function [sampleNames, keepGoodData, exportCurves] = considerMAVENv14(CSVfileFromMAVEN,sampleInfoFile,setQuality,pruneData,ionMode_SRMfile)
%1: CSVfileFromMAVEN: the name of the CSV file from MAVEN
%2: sampleInfoFile: the details on where the file names for the samples, 
%blanks, and standards 
%3: setQuality: the required level for a peak in the standard curve to be
%considered as good. Have been using 0.5, but might want to be more lenient
%at times
%working on using the detailed information from MAVEN to construct
%calibration curves. Want to set specific criteria before I will allow a
%peak to be used in a calibration curve. These criteria will be: (1) the
%quality score from MAVEN, (2) the signal:noise ratio of the peak in MAVEN,
%and (3) something that considers the range of values found in the samples
%4. pruneData, set to 1 if I want to only consider files marked as goodData
%5. ionMode is positive/negative/both --> three runs on the TSQ
%6. SRMfile is the name of the CSV file with the SRM list, one per ion mode
%in the sample info file
%KL 4/21/2014; KL 5/19/2014 (v2, added check for confirm ion; 
%KL 6/10/2014 (v3, need to treat d2-biotin pos/neg separately bc don't want
%standard curve).
%KL 7/28/2014 (v5 adds the name of the sample information file to the input
%rather than relying on changing it in the middle of this function)
%KL 10/22/2014, v6 is updated to use a new list of SRMs and the extended
%calibration curve range (up to 1000 ng/ml now)
%KL 5/19/2016 v7 modify for Alm samples, but making other changes to error
%checking
%KL 6/29/2016 v8 broadening the error correction after the Alm samples
%then realized I need a special case bc standard cruve in odd order
%KL 10/20/2016 for the new UPLC method need some different files for positive
%and negative ion mode, also adding in Winn's QC for confirm ions
%KL 2/16/2017 for Chisholm/Pro PLimited1 samples, allow poorer quality data
%for methyl-pentanoic acids
%KL 1/16/2018 remove spike from the list of samples to be used when
%considering the range of concentrations to use in the standard curve. Too
%often the spike is well above the rest of the samples
%KL 2/21/2020 use SRM definitions from riMAVEN12
%KL 3/12/2021 updating to use cleaner names from Melissa - should minimize
%special cases. Also need to add in Mix8 with both positive and negative
%ion mode data in one run; also add ability to keep points for standard
%curves bc need that information
%%get started:
%%%%use the SRM file to setup the dataset array to be filed with data, only
%%%%using the main compounds and not the confirm ones

% warning('off','MATLAB:table:ModifiedAndSavedVarnames')
% warning('off','stats:dataset:ModifiedVarnames')
warning('off')

switch ionMode
    case 'positive'
        %maybe this will be easier to read:
        special = table();
        special.cases(1) = {'2-deoxyinosine_Na pos'}; special.match(1) = {''};
        special.cases(2) = {'asparagine'}; special.match(2) = {'asparagine confirm 2'};
        special.cases(3) = {'isoleucine'}; special.match(3) = {'isoleucine confirm 2'};
        
    case 'negative'
        special = table();
        special.cases(1) = {'desthiobiotin neg'}; special.match(1) = {''};
        special.cases(2) = {'3-methyl-2-oxobutanoic acid'}; special.match(2) = {''};
        special.cases(3) = {'4-methyl-2-oxopentanoic acid'}; special.match(3) = {''};
        special.cases(4) = {'3-3-dimethyl-2oxobutanoic acid'}; special.match(4) = {''};
        special.cases(5) = {'2-methyl-4-oxopentanoic acid'}; special.match(5) = {''};
        special.cases(6) = {'isethionic acid '}; special.match(6) = {'isethionic acid confirm 1'}; 
        special.cases(7) = {'2oxohexanoic acid'}; special.match(7) = {''};
        
    case 'both'
        special = table();
        special.cases(1) = {'indole-3-acetic acid d7'}; special.match(1) = {'indole-3-acetic acid d7confirm'};
        special.cases(2) = {'3-methyl-2-oxopentanoic acid'}; special.match(2) = {''};
        
end

srm = dataset('File',SRMfile,'delimiter',',');

%easiest to find the confirm compounds
r = regexp(srm.compound,'confirm');
gm = cellfun(@isempty,r); %will be 1 when the compound is NOT labeled as confirm
k = find(gm==1);

compoundList = mat2dataset(srm.compound(k,1),'VarNames',{'name'}); 
compoundList.indexMain = 0;
compoundList.indexConfirm = 0;
clear r gm k

%go through each compound, find the index for the main peak and the confirm 
for a = 1:length(compoundList);
    %change syntax 3/16/2021
    %sn = strcmp(compoundList.name(a),srm.compound);
    %ks = find(sn==1);
    ks = strmatch(compoundList.name(a),srm.compound,'exact');
    compoundList.indexMain(a,1) = ks;
    
    h = strmatch(compoundList.name(a),special.cases);

    if isempty(h)
        %these will match as is with no need to dig into the special cases
        %3/12/2021 while I am working on this, just make this broader - match
        %anything with the name and confirm ...
        %sn = contains(srm.compound,'confirm') & contains(srm.compound,compoundList.name(a));
        tn = strcat(compoundList.name(a),' confirm');
        sn = strcmp(tn,srm.compound);
        ks = find(sn==1);
%         try 
            compoundList.indexConfirm(a,1) = ks;
%         catch
%             fprintf('here')
%         end
        clear ks tn sn
    elseif ~isempty(h) && ~isempty(special.match{h})
        %setup the name to use in special - go find that
        s2 = strcmp(special.match(h),srm.compound);
        ks2 = find(s2==1);
        compoundList.indexConfirm(a,1) = ks2;
        clear s ks s2 ks2
    elseif ~isempty(h) && isempty(special.match{h})
        %no confirm ion
        compoundList.indexConfirm(a,1) = NaN;
    end
    clear sn ks h
    
end
clear a


% %%%%% import the CSV file from MAVEN into a dataset array
% %%%%% import the CSV file from MAVEN into a dataset array

%now go ahead and use the dataset function to open the file
importedPeakList = dataset('File',CSVfileFromMAVEN,'delimiter',',');

%clean up
clear CSVfileFromMAVEN tempFile fid tline fidOut tline ans

% warning('off', 'stats:dataset:genvalidnames:ModifiedVarnames');

%%change this if the sampleInfoFile is a CSV file:
sampleInfo = dataset('XLSFile',sampleInfoFile);
% sampleInfo = dataset('File',sampleInfoFile,'delimiter',',','HeaderLines',1) %works for CSV files
clear CSVfileFromMAVEN

%can have cases where I have decided files in the TSQ list should not be
%processed
if pruneData
    k = find(sampleInfo.goodData==0);
    sampleInfo(k,:)=[];
    clear k
end

%setup empty matrix for the data and a separate matrix for error
%how many possible standards are there, and where will they be?
%KL 10/20/2016 now need to look for positive or negative set
%KL 3/15/2021 updated to pull concentrations from file
switch ionMode
    case 'negative'
        kStandard = find((strcmp(sampleInfo.SampleType,'Std Bracket')==1) & (strcmp(sampleInfo.ionMode,'negative')==1));           
        %setStandardConcentrations = [0.5 1 5 10 25 50 100 250 500 1000]';
    case 'positive'
        kStandard = find((strcmp(sampleInfo.SampleType,'Std Bracket')==1) & (strcmp(sampleInfo.ionMode,'positive')==1));   
        %setStandardConcentrations = [0.5 1 5 10 25 50 100 250 500 1000]';
    case 'both'
        kStandard = find((strcmp(sampleInfo.SampleType,'Std Bracket')==1) & (strcmp(sampleInfo.ionMode,'both')==1));          
end

setStandardConcentrations = sampleInfo.concentration_ngML(kStandard); %remember - dataset...not table
        
standardNames = sampleInfo.FileName(kStandard);
nStandards = length(kStandard);

kSample = find(strcmp(sampleInfo.SampleType,'Unknown')==1);
sampleNames = sampleInfo.FileName(kSample);
sampleNames_KLworking = sampleInfo.SampleType(kSample);
nSamples =length(kSample);
clear kStandard kSample

goodData(1:length(compoundList),nSamples) = NaN; 
goodDataError = goodData;

warning('off', 'stats:dataset:subsasgn:DefaultValuesAddedVariable');
compoundList.r2_line = 0;
compoundList.slope = 0;
compoundList.intercept = 0;
compoundList.SDslope = 0;
compoundList.SDintercept = 0;
compoundList.nPoints = 0;


%Liz wants 5 points in the standard curve (a/o 4/24/2014), so set that here
nRequired = 5;

%go through one compound at a time, (1) make the standard curve, (2) use that to
%calculate the areas for each compound for each sample, (3) then go find the
%confirm ion for each compound in each sample and make sure I like where
%things are

%note that I want to skip the d2-Biotin pos/neg because that does not have
%a standard curve (though there is some d2Biotin in the higher levels of
%the standards.) However, I just need the raw values for every sample
biotin = {'d2-biotin neg','d2-biotin pos'};
exportCurves = table(); %adding 3/24/2021 to send this info outside fxn

for a = 1:length(compoundList);
       
    %keep the points for the standard curve - need this outside the function
    exportCurves.mtabName(a) = compoundList.name(a);
    exportCurves.concentration{a} =  setStandardConcentrations;

    if strcmp(compoundList.name(a),'trehalose_387 neg')
        %stop here for troubleshooting as needed
        fprintf('here')
        %compoundList.name(a)
    end
    
    clear xdata ydata %here to make checking out one compound easier
    k = strmatch(compoundList.name(a),importedPeakList.compoundId,'exact');
           
    if ~isempty(k)
        %can have cases where nothing good was found. If k is NOT empty,
        %found good data

        smallDS = importedPeakList(k,:); clear k
        [c ia ib] =intersect(smallDS.sample,sampleInfo.FileName);
        smallDS.sampleType(ia,1) = sampleInfo.SampleType(ib,1);
        clear c ia ib  

        [~, idxDS, idxStandards] = intersect(smallDS.sample,standardNames);

        %what is the average value in the blanks? Need this for two reasons,
        %(1) to get a zero value for the standard curve and (2) to see if the
        %values in the samples are more/less than what is in the blanks
        kb = find(strcmp(smallDS.sampleType,'Blank')==1);
        %for now, using AreaTop, might play around with that later
        meanBlank = mean(smallDS.peakArea(kb)); clear kb

        xdata = setStandardConcentrations;
        ydata(1:length(xdata),1) = NaN;
        quality = ydata;
        %get all possible values from the standard curve
        ydata(idxStandards) = smallDS.peakArea(idxDS);
        quality(idxStandards) = smallDS.quality(idxDS);

        %original = [xdata ydata quality];
        exportCurves.show{a} = [xdata ydata quality];
        
        % now prune the standard data to only include good data:
        %let's just set the bad values to NaN here so that I don't have to keep
        %resizing everything...so, find the BAD data
        
        %for some compounds, allow lower quality scores...new acids in
        %particular are not so great; KL 4/9/2021
        lowQuality = {'4-methyl-2-oxopentanoic acid',...
            '3-methyl-2-oxopentanoic acid',...
            '2-methyl-4-oxopentanoic acid',...
            '3-3-dimethyl-2oxobutanoic acid',...
            'sucrose_341 neg',...
            'trehalose_387 neg'}';
        if isequal(sum(strcmp(compoundList.name(a),lowQuality)),1)
            kBad = find(quality < 0.05);
            dBad = find(smallDS.quality < 0.05);
        else 
            kBad = find(quality < setQuality);
            dBad = find(smallDS.quality < setQuality);
        end         
            
        ydata(kBad) = NaN;
        xdata(kBad) = NaN; 
        clear quality
        
        %from Winn (10/2016): remove low quality peaks from the data        
        smallDS.peakArea(dBad) = NaN;

        %put the blank at the beginning...remember that this will treat all
        %blanks equally...this may be a bad thing
        xdata = cat(1,0,xdata); 
        
        ydata = cat(1,meanBlank,ydata); 
        %clear meanBlank
        clear kBad idxDS idxStandards      

        %%%before I go ahead and calculate the values for each sample, make
        %%%sure that the confirm ion is present at the same RT as the main ion.
        %%%adding this in following Melissa's suggestion, 5/19/2014
        i = compoundList.indexConfirm(a);
        %if there is no confirm ion, can't do the check...desthiobiotin neg
        %has no confirm. For now, allow that to pass as if it is always OK
        if ~isnan(i)
            k = strmatch(srm.compound(i),importedPeakList.compoundId,'exact');

            %k will be empty if no confirm ions were found  
            if ~isempty(k) 
                smallDSconfirm = importedPeakList(k,:); clear k
                smallDS.cfRT(1:size(smallDS,1),1) = NaN;

                [c ia ib] = intersect(smallDSconfirm.sample,smallDS.sample);
                %now compare the retention times for the main ion and the confirm
                smallDS.cfRT(ib,1) = smallDSconfirm.rt(ia);  
                %from Winn (10/2016), require the confirm to also meet a 
                % QC check (for now I am using setQuality/10):
                %Winn required 0.1, but with the UPLC data this cuts most
                %thymidine confirm ions
                %also organize the confirm ion QC data
                smallDS.cfquality(ib,1) = smallDSconfirm.quality(ia);
                smallDS.cfsig(ib,1) = smallDSconfirm.signalBaseLineRatio(ia);
                %set the requirements for the quality of the confirm ion
                %peak and convert the confirm RT to NaN for those that fall
                %below the threshold:
                lowqc = smallDS.cfquality < (setQuality/10) | (smallDS.cfsig <= 1);
                smallDS.cfRT(lowqc) = NaN;              
                d = abs(minus(smallDS.cfRT,smallDS.rt));

                %require that the confirm and the main peak have retention times within
                %10 seconds (remember to convert to minutes!)
                %reqRT = 10/60; %this was what I used for the ventDOM project
                reqRT = 12/60; %MCKS suggested 5/22/2014 that this would be better

                k2 = find(d > reqRT | isnan(d)); %will be a NaN if peak is not found

                %for the peaks that fail the confirm retention time check, set them
                %equal to NaN
                smallDS.peakArea(k2) = NaN;
                clear k2 d c ia ib
            elseif isempty(k)
                %no, will deem this as not good data bc no confirm ions were
                %found
                smallDS.peakArea(:) = NaN;
            end
        end

        %remember, will also have cases where no data were found for select samples
        %so need to setup the spacers in there are well
        [c ia ib] = intersect(smallDS.sample,sampleNames);
        tData(1:length(sampleNames),1) = NaN;
        tData(ib) = smallDS.peakArea(ia); 
        clear c ia ib
        full=tData;
        
        % 1/16/2018 adding in the ability to ignore the pooledSpike when considering
        % what to include in the standard curve (slight hack)
        % need the data that is not from the pooled spike
        su = strcmp(sampleInfo.sampleDetails,'Unknown');
        ksu = find(su==1);
        [c ia ib] = intersect(sampleInfo.FileName(ksu),smallDS.sample);
        %6/26/2018 can have the case where no unknowns get out of MAVEN...
        if ~isempty(c)
            tData_unknownsOnly = smallDS.peakArea(ib);
        else
            tData_unknownsOnly = NaN;
        end
        clear su ksu c ia ib
                    
        %now go and start with the special case for the d2-Biotin standards
        if sum(strcmp(compoundList.name(a),biotin)>0)
            %this is one of the biotin standards, treat separately
            %just need the peak areas, not converted to anything
            goodData(a,:) = tData;
            goodDataError(a,:) = NaN;

        else %not biotin, move on
            %try adding this...if the value is less than the meanBlank, change
            %that to NaN...adding 5/18/2016
            k = find(tData<meanBlank);
            tData(k) =NaN;
            clear k meanBlank

            %another option is to not allow values below the smallest 'good'
            %value in the measured peak areas...add this 5/18/2016
            k = find(tData< min(ydata(2:end)));
            tData(k) = NaN;
            clear k
        
            %need to deal with the idea of how big to allow the curve to be
            %and, what is the max value in my samples? should probably have at least one point above that
            m = max(tData_unknownsOnly); 
            %if all the unknowns fail the quality check, this next step
            %will fail. Haven't seen this until now (6/26/2018)
            if isnan(m)
                %easiest to make kMax empty
                kMax = [];
            else
                kMax = find(ydata <= m);
            end            
            clear tData_unknownsOnly

            if ~isempty(kMax)
                %have at least one point on the curve
                if isequal(kMax(end),length(ydata)) 
                    %already at the end of the standard curve...so use all the points
                    %do nothing...but send up a flag since the data are above
                    %the standard curve
                    disp([compoundList.name(a) ' is above the standard curve'])
                    %fprintf('here')
                elseif isequal(kMax(end)+1,length(ydata));
                    %only one more above the points in the standard curve, use all the
                    %points
                elseif isequal(kMax,1);
                    %data are at the low end of the standard curve, but let's require 
                    %more points above my data to get a reasonable curve...
                    xdata = xdata(1:nRequired);
                    ydata = ydata(1:nRequired);
                elseif length(kMax)+2  < nRequired
                    % use the number of points sent in nRequired
                    ydata = ydata(1:nRequired);
                    xdata = xdata(1:nRequired);
                elseif length(kMax) + 1 < nRequired
                    %use the standard curve to one point beyond the range of my
                    %samples
                    ydata = ydata(1:kMax(end)+1);
                    xdata = xdata(1:kMax(end)+1);
                else
                    %use the standard curve to one point beyond the range of my
                    %samples
                    ydata = ydata(1:kMax(end)+1);
                    xdata = xdata(1:kMax(end)+1);

                end
            elseif isempty(kMax)
                %all of the points in the standard curve are higher than what was
                %measured in the samples
                ydata = ydata(1:nRequired);
                xdata = xdata(1:nRequired);
            end
            clear kMax m

            %need at least three points to make a curve AND get the error estimates
            try
                show = [xdata ydata];
            catch 
                error('Something wrong, xdata and ydata should be the same size')
            end
            i = isnan(show);
            sfmi = sum(i,2);
            k = find(sfmi==0);
            xdata = xdata(k);
            ydata = ydata(k); 
            clear show i sfmi k
            %this will be helpful bc will show where I had <2 points  
            %remember that this also takes into account the rules I set above about
            %how wide to make the standard curve
            compoundList.nPoints(a) = length(ydata); 
                    
        if length(xdata)>2
                dataOut = getErrors(xdata,ydata); %errors for the standard curve
                [calcError, calcConc] = useErrors(dataOut,tData); %then calculate the concentrations

                %add in a check, if the slope is negative, this is garbage
                if dataOut.slope > 0;          
                    %this will be the same number of rows as unCompounds
                    %the number of columns will match the number of unknown samples
                    goodData(a,:) = calcConc; 
                    goodDataError(a,:) = calcError; %can get percent by calcError./calcConc

                    compoundList.slope(a) = dataOut.slope;
                    compoundList.intercept(a) = dataOut.intercept;
                    compoundList.SDslope(a) = dataOut.SDslope;
                    compoundList.SDintercept(a) = dataOut.SDintercept;
                    compoundList.r2_line(a) = dataOut.r2;
                else
                    goodData(a,:) = NaN;
                    goodDataError(a,:) = NaN;

                    compoundList.slope(a) = NaN;
                    compoundList.intercept(a) = NaN;
                    compoundList.SDslope(a) = NaN;
                    compoundList.SDintercept(a) = NaN;
                    compoundList.r2_line(a) = NaN;
                end

            else
                %not enough points to make a standard curve
                goodData(a,:) = NaN;
                goodDataError(a,:) = NaN;

                compoundList.slope(a) = NaN;
                compoundList.intercept(a) = NaN;
                compoundList.SDslope(a) = NaN;
                compoundList.SDintercept(a) = NaN;
                compoundList.r2_line(a) = NaN;
            end
             %%%DE-BUGGING HERE
             %compoundList.name(a)
             %put breakpoint at the next line and uncomment out the
             %compoundList.name(a) line above if troubleshooting one
             %compound at a time (5/19/2016)
             clear dataOut calcError calcConc tData

             clear xdata ydata smallDS
        end %this is the end of the IF biotin 
 
    end
end %end of calculating the curves and data
clear a

ds2 = dataset(goodData);
ds3 = dataset(goodDataError);
keepingAll = cat(2,compoundList,ds2,ds3);
clear ds1 ds2 ds3 compoundList goodData goodDataError

%remove some unneeded variables:
keepingAll.indexMain = [];
keepingAll.indexConfirm = [];

%perhaps do a little pruning to provide a dataset array with only the data
%that are good from the criteria above and are not overlapping with
%zero...
i = isnan(keepingAll.SDintercept);
k = find(i~=1);

keepGoodData = keepingAll(k,:);
%here, we need to consider on a sample by sample basis and not make
%decisions based on the entire set for each compound

for a = 1:size(keepGoodData,1);
    %check that this is not biotin...
    if sum(strcmp(keepGoodData.name(a),biotin)>0)
        %keep all the biotin data, so can do nothing here
    else
        for aa = 1:size(keepGoodData.goodData,2);
            tD = keepGoodData.goodData(a,aa);
            tE = keepGoodData.goodDataError(a,aa);

            %can have a few options.
            if tD < 0 %easiest: sample is less than zero
                keepGoodData.goodData(a,aa)=0;
                keepGoodData.goodDataError(a,aa) = 0;
            elseif tD - tE < 0 %does the error window include zero?
                %what is the window around it,
                keepGoodData.goodData(a,aa)=0;
                keepGoodData.goodDataError(a,aa) = 0;
%             elseif tE./tD*100 > 66 %is the error percent above 66%?
%                 %added 5/18/2016 bc getting too many things that have 
%                 %values that get calculated, but the error is high cfd
%                 %to the measured value
%                 keepGoodData.goodData(a,aa)=0;
%                 keepGoodData.goodDataError(a,aa) = 0;

            end
            clear tD tE
        end
        clear aa
    end
end
clear a

%can also have the case where I am now left with zeros and NaNs only...
i = isnan(keepGoodData.goodData);
for a = 1:size(i,1);
    td = keepGoodData.goodData(a,:);
    ki = find(i(a,:)==1);
    k = find(i(a,:)~=1);
    ts = sum(td(k));
    if ts==0
        keepGoodData.goodData(a,ki) = 0;
    end
    clear td ki k ts
end
clear a

%now go ahead and delete the rows where all the datapoints are zero...no
%sense making my life more difficult
fm = logical(keepGoodData.goodData~=0);
sfmc = sum(fm,2);
k = find(sfmc==0);

keepGoodData([k],:) = [];
    
clear fm sfmc k

%%%%put the internal functions here at the end
%%%%put the internal functions here at the end

    function dataOut = getErrors(xdata,ydata);
    %function dataOut = getErrors(xdata,ydata);
    %From this web site:
    %http://terpconnect.umd.edu/~toh/spectrum/LeastSquaresMatlab.txt
    %KL modifying 4/21/2014

    x = xdata;
    y = ydata;
    % Simple Matlab script for calculating the first-order least-square fit of y vs x,
    % including the Slope and Intercept and the predicted standard deviation of
    % the slope (SDSlope) and intercept (SDIntercept).

    NumPoints=length(x);
    Sxx = sum((x-mean(x)).^2);
    Syy = sum((y-mean(y)).^2);
    Sxy = sum((x-mean(x)).*(y-mean(y)));
    Slope = Sxy./Sxx;
    Intercept = mean(y)-Slope*mean(x);

    Sy = sqrt((Syy-Slope^2*Sxx)/(NumPoints-2));

    SDslope = Sy/sqrt(Sxx);
    SDintercept = Sy*sqrt(1./(NumPoints-(sum(x).^2)./sum(x.^2)));

    r2 = 1 - ((Syy-Slope^2*Sxx) ./Syy);

    %data to send out of this function (when it is a function)
    dataOut.slope = Slope;
    dataOut.intercept = Intercept;
    dataOut.SDslope = SDslope;
    dataOut.PercentSlopeError = SDslope./Slope;
    dataOut.SDintercept = SDintercept;
    dataOut.PercentInterceptError = SDintercept./Intercept;
    dataOut.r2 = r2;

    end %end of getErrors as a function


    function [calcError, calcConc] = useErrors(myErrorData,measuredSample)
    %function [calcError, calcConc] = useErrors(Slope,Intercept,SDslope,SDintercept,measuredSample)
    %use the errors on the line to get the errors on the samples actually
    %measured
    %KL 4/21/2014
    Intercept = myErrorData.intercept;
    Slope = myErrorData.slope;
    SDintercept = myErrorData.SDintercept;
    SDslope = myErrorData.SDslope;

    %calculated concentrations from my hypothetical list
    calcConc = (measuredSample - Intercept)./Slope;

    %apply to my list of hypothetical unknowns, split this up to make
    %it easier to keep track of where the parentheses etc. are
    fSQ = (SDintercept./(measuredSample - Intercept)).^2 + (SDslope./Slope).^2;
    calcError = calcConc .* sqrt(fSQ);
    errorPercent = calcError./calcConc*100;

    end %end of useErrors as a function

end %end of considerMAVEN as a function



            


