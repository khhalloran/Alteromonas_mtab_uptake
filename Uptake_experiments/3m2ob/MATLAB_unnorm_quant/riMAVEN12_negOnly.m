%get the MAVEN data and the sample information squared away
%this code will work for UPLC data where two runs: pos and neg mode
%allow poor quality scores for the new ''methyl-*tanoic' acids
%KL 10/20/2016; KL 11/14/2017 add in the ability to merge based on names
%KL 1/20/2018 - around line 245 is code that works here (MATLAB 2013b, but
%will not work in MATLAB 2017a). 
%KL 6/4/2019 - doing some housecleaning
clear

%%%set the file names up front:
%MAT file for MATLAB
NameOfFile = 'KHU7_TSQ.2021.05.24.mat';

%setup some variables for the TSQ processing:
setQuality = 0.2; %quality score comes from MAVEN analysis
warning('off', 'MATLAB:codetools:ModifiedVarnames');

%where is sequence file from the TSQ?
wDir = '/Users/kathrynhalloran/Documents/MATLAB/Uptake experiments/3m2ob_7_v2/sequence_fromMethods';
fName = 'mtab_KHU7_2021_0419.xls';
sampleInfoFile = [wDir filesep fName];
clear wDir

%where are the SRM lists that will be needed for MATLAB (same files as used
%in MAVEN, but modified from the versions exported out of the TSQ)?
sDir = '/Users/kathrynhalloran/Documents/MATLAB/Uptake experiments/3m2ob_7_v2/SRM_list';
%srmFile_pos = [sDir filesep 'Pos_srms_uplc_forMAVEN.2019.04.26.csv'];
srmFile_neg = [sDir filesep 'Neg_UPLC_2021.05.06_forMAVEN_KHU7.csv']; 
clear sDir

%Since there is no comma at the end of the headerline, MATLAB will crash
%with the import. This means we have to use some low level input/output
%commands to read in the text files
%...read in the first line, make sure it ends in a comma, 
%close the file and then open it again...

%start by setting up the list of CSV files exported out of MAVEN:
% readFile1_pos = 'TSQ_Weber_coral_pos_part1.2019.05.02.csv';
% readFile2_pos = 'TSQ_Weber_coral_pos_part2.2019.05.02.csv';
% readFile3_pos = 'BIOSSCOPE_TSQexport_pos_part3.2018.01.04.csv';

readFile1_neg = 'KHU7_peaklist_v2.csv';
% readFile2_neg = 'TSQ_Weber_coral_neg_part2.2019.04.29.csv';
% readFile3_neg = 'TSQ_Weber_coral_neg_part3.2019.04.29.csv';
% readFile4_neg = 'TSQ_Weber_coral_neg_part4.2019.05.02.csv';
% readFile5_neg = 'TSQ_Weber_coral_neg_part5.2019.04.29.csv';


%%now do the negative ion mode files
%%now do the negative ion mode files
tempFile = 'temp.csv';
fid1 = fopen(readFile1_neg);
tline = fgetl(fid1);
if ~strcmp(tline(end),',')
    tline = strcat(tline,','); %put the comma at the end of the headerline
end
%open up the temp file, and set it to be write-able
fidOut = fopen(tempFile,'w');
fprintf(fidOut, '%s\n',tline);

tline = fgetl(fid1); %get the next line of the file
while ischar(tline)
    fprintf(fidOut, '%s\n', tline);
    tline = fgetl(fid1);    
end
fclose(fid1);

%now get the second file...since the headerline will be the same can skip it
% fid2 = fopen(readFile2_neg);
% tline = fgetl(fid2);
% tline = fgetl(fid2); %get the next line of the file
% while ischar(tline)
%     fprintf(fidOut, '%s\n', tline);
%     tline = fgetl(fid2);    
% end
% fclose(fid2);
% 
% if 1
%     %now get the third file...since the headerline will be the same can skip it
%     fid3 = fopen(readFile3_neg);
%     tline = fgetl(fid3);
%     tline = fgetl(fid3); %get the next line of the file
%     while ischar(tline)
%         fprintf(fidOut, '%s\n', tline);
%         tline = fgetl(fid3);    
%     end
%     fclose(fid3);
%     
%     %now get the fourth file...since the headerline will be the same can skip it
%     fid4 = fopen(readFile4_neg);
%     tline = fgetl(fid4);
%     tline = fgetl(fid4); %get the next line of the file
%     while ischar(tline)
%         fprintf(fidOut, '%s\n', tline);
%         tline = fgetl(fid4);    
%     end
%     fclose(fid4);
%     
%     
%     %now get the fifth file...since the headerline will be the same can skip it
%     fid5 = fopen(readFile5_neg);
%     tline = fgetl(fid5);
%     tline = fgetl(fid5); %get the next line of the file
%     while ischar(tline)
%         fprintf(fidOut, '%s\n', tline);
%         tline = fgetl(fid5);    
%     end
%     fclose(fid5);
% end


fclose(fidOut);
clear fid* tline ans readFile*_neg 

%%%now move on with the TSQ data processing
[neg.sNames neg.kgd] = considerMAVENv12(tempFile,sampleInfoFile,srmFile_neg,setQuality,'negative');
clear fName tempFile srmFile_neg
clear setQuality setErrorP


%%now merge the pos and neg mode data from the two different UPLC runs
%%now merge the pos and neg mode data from the two different UPLC runs

%get all the metabolite names...should be a unique list even when I merge
%across positive and negative ion mode, but check that
mtabNames = neg.kgd.name;
% mtabNames = sort(cat(1,neg.kgd.name,pos.kgd.name));
% if length(unique(mtabNames)) ~= length(mtabNames)
%     error('Something is wrong - duplicate names in the list of metabolites')
% end

%for the pooled samples (and perhapds others), I will have duplicate sets 
%of names with either _pos or _neg appended; 
tInfo = readtable(sampleInfoFile);
clear sampleInfoFile

% %%first, go through and iterate through the pooled samples
% %%to provide numbers for these (otherwise will have duplicate
% %%names)
% %%NOTE: update names of pooled samples here
% s = strcmp(tInfo.SampleName,'Coral_Pool_pos');
% ks = find(s==1);
% for a = 1:length(ks);
%     t = tInfo.SampleName(ks(a));
%     tInfo.SampleName(ks(a)) = strcat('p',num2str(a),t);
%     clear t
% end
% clear a ks a
    
%%first, go through and iterate through the pooled samples
%%to provide numbers for these (otherwise will have duplicate
%%names)
%%NOTE: update names of pooled samples here
s = strcmp(tInfo.SampleName,'KHU7_Pool');
ks = find(s==1);
for a = 1:length(ks);
    t = tInfo.SampleName(ks(a));
    tInfo.SampleName(ks(a)) = strcat('p',num2str(a),t);
    clear t
end
clear a ks a

%before I dive into the unknowns, remove anything that has goodData = 0
k = find(tInfo.goodData==0);
tInfo(k,:) = [];
clear k

% %now find the Unknown...should have the same number for positive and
% %negative ion mode bc have pruned out the different QC samples already
% s = strcmp(tInfo.SampleType,'Unknown');
% sp = strcmp(tInfo.ionMode,'pos');
% ksp = (find(s==1 & sp==1));
% sn = strcmp(tInfo.ionMode,'neg');
% ksn = (find(s==1 & sn==1));
% 
% if ~isequal(length(ksp),length(ksn))
%     error('Something wrong, these should be the same length')
% end
% clear s sp sn ksp ksn

%%parse out the names. Use this to figure out the unique samples and setup
%%a new matrix that I can propagate with the metabolites from both positive
%%and negative ion mode. Bit of a hack, and growing worse.
nrow = size(tInfo,1);
tInfo.type = repmat({''},nrow,1);
tInfo.cName = repmat({''},nrow,1);
%examples of additional columns used in the BIOS-SCOPE project
% tInfo.cruise = repmat({''},nrow,1);
% tInfo.cast = zeros(nrow,1);
% tInfo.niskin = zeros(nrow,1);
% tInfo.depth = zeros(nrow,1);
% tInfo.addedInfo = repmat({'none'},nrow,1);

for a = 1:nrow;
    if strcmp(tInfo.SampleType{a},'Unknown') %only do unknowns      
        one = tInfo.SampleName{a};
        r_pooled = regexp(one,'KHU7_Pool');
        r_spiked = regexp(one,'spike');

            if r_spiked
                %pooled with 500 ng/ml spike
                if 0
                    %keep it
                    tInfo.type(a) = {'spiked'};
                    tInfo.cName(a) = {'spiked'};
                else
                    %skip
                end
            elseif r_pooled
                %pooled sample
                tInfo.type(a) = {'pooled'};
                %put the number of this pooled sample into 'addedInfo'
                r_nL = regexp(one,'p'); %lower case
                r_nU = regexp(one,'P'); %upper case
                %tInfo.addedInfo(a) = {one(r_nL+1 : r_nU-1)};
                tInfo.addedInfo(a) = {'pooled'};
                tInfo.cName(a) = {one(1:r_nU-1)};
            else
                %actual sample
                tInfo.addedInfo(a) = {'sample'}; %redundant...
                tInfo.cName(a) = {one(1:end-4)};
                %fprintf('here')
            end
        clear one r_* under
    end
end
clear a nrow

sInfo = table;
sInfo.cName = unique(tInfo.cName);
%the first row of this will be empty, delete that
if isequal(sInfo.cName(1),{''});
    sInfo(1,:) = [];
end


%now make an empty matrix for the data...will be all numbers so no need for
%special format
mtabData = zeros(size(mtabNames,1),size(sInfo,1));
% %need to track some additional details:
% mtabDetails = table();

%get the index for rows for positive AND negative mtabs:
%[c idx_posNew idx_posOld] = intersect(mtabNames,pos.kgd.name);
[c idx_negNew idx_negOld] = intersect(mtabNames,neg.kgd.name);

% mtabDetails.mode(idx_posNew,1) = {'pos'};
% mtabDetails.mode(idx_negNew,1) = {'neg'};

% sInfo.runOrder_pos(:,1) = 0;
sInfo.runOrder_neg(:,1) = 0;

% sInfo.FileName_pos(:,1) = {''};
sInfo.FileName_neg(:,1) = {''};

% for a = 1:size(sInfo,1);
%     s = strcmp(sInfo.cName(a),tInfo.cName);
%     ks = find(s==1);
%     if length(ks) ~= 2
%         error('Something is wrong, should be two of each')
%     end
%     
%     %some variant of this:
%     for aa = 1:2
%         %propagate sInfo with the cast/depth/etc. information, only do once
% %         if aa == 1
% %             sInfo.type(a) = tInfo.type(ks(aa));
% %             sInfo.cName(a) = tInfo.cName(ks(aa));
% %             sInfo.cruise(a) = tInfo.cruise(ks(aa));
% %             sInfo.cast(a) = tInfo.cast(ks(aa));
% %             sInfo.niskin(a) = tInfo.niskin(ks(aa));
% %             sInfo.depth(a) = tInfo.depth(ks(aa));
% %             sInfo.addedInfo(a) = tInfo.addedInfo(ks(aa));
% %         end
%         
%         im = tInfo.ionMode{ks(aa)};
%         if isequal(im,'pos')
%             tName = tInfo.FileName(ks(aa));
%             sInfo.FileName_pos(a,1) = tName;
% 
%             [c ia tIdx] =intersect(tName,pos.sNames);
%             mtabData(idx_posNew,a) = pos.kgd.goodData(idx_posOld,tIdx);
%             clear c ia tIdx tName
%             
%         elseif isequal(im,'neg')
%             tName = tInfo.FileName(ks(aa));
%             sInfo.FileName_neg(a,1) = tName;
% 
%             [c ia tIdx] =intersect(tName,neg.sNames);
%             mtabData(idx_negNew,a) = neg.kgd.goodData(idx_negOld,tIdx);
%             clear c ia tIdx tName
%         else 
%             error('Something wrong')
%         end
%         clear im
%     end
%     clear aa s ks        
% end
% clear a
% 
% clear idx_* tInfo
% 
% % for a = 1: size(sInfo,1)
% %     %do positive ion mode first
% %     gc = sInfo{a,'FileName_pos'}{:}; %added {:} to deal with table output
% %     t = regexp(gc,'_');
% %     if ~isempty(t)
% %         sInfo.runOrder_pos(a,1) = str2num(gc(t(end)+1:end));
% %     else
% %         sInfo.runOrder_pos(a,1) = NaN;
% %     end
% %     clear gc t
% %     
% %     %then negativeion mode first
% %     gc = sInfo{a,'FileName_neg'}{:}; %added {:} to deal with table output
% %     t = regexp(gc,'_');
% %     if ~isempty(t)
% %         sInfo.runOrder_neg(a,1) = str2num(gc(t(end)+1:end));
% %     else
% %         sInfo.runOrder_neg(a,1) = NaN;
% %     end
% %     clear gc t
% % end
% % clear a
%  
% %making decisions about pos/neg mtabs; this is done either based on past
% %experience with compounds, decisions made within a project, by plotting
% %the positive and negative ion mode data together, or looking at the plots
% %in MAVEN. The first time this code is run, I allow everything to go
% %through and then I later rerun the code to make decisions to leave one
% %version for each compound
% if 0
%     toDelete = {'5deoxyadenosine neg','NAD neg','adenosine 5''-monophosphate neg',...
%         'adenine neg','biotin neg','cytidine neg','desthiobiotin neg',...
%         'folic acid neg', 'guanosine neg','n-acetyl muramic acid neg',...
%         'pantothenic acid pos','s-(5''-adenosyl)-L-homocysteine neg',...
%         'sn-glycerol 3-phosphate pos ','xanthine neg','xanthosine pos','hemin II',...
%         '2deoxyguanosine neg','2deoxycytidine neg','3deoxyguanosine neg',...
%         'aspartic acid neg'};
% 
%     [c ia ib] = intersect(toDelete,mtabNames);
% 
%     mtabNames(ib,:)=[];
%     mtabDetails(ib,:)=[];
%     mtabData(ib,:)=[];
%     clear c ia ib toDelete
% end
% 
% %manually house clean a few compounds:
% s = strcmp(mtabNames,{'D-Ribose 5-phosphate'});
% ks = find(s==1);
% mtabData(ks,:) = [];
% mtabDetails(ks,:) = [];
% mtabNames(ks,:) = [];
% clear s ks
% 
% s = strcmp(mtabNames,{'D(-)3-phosphoglyceric acid'});
% ks = find(s==1);
% mtabData(ks,:) = [];
% mtabDetails(ks,:) = [];
% mtabNames(ks,:) = [];
% clear s ks
% 
% s = strcmp(mtabNames,{'glutamic acid d3 pos'});
% ks = find(s==1);
% mtabData(ks,:) = [];
% mtabDetails(ks,:) = [];
% mtabNames(ks,:) = [];
% clear s ks
% 
% s = strcmp(mtabNames,{'phenylalanine 13C'});
% ks = find(s==1);
% mtabData(ks,:) = [];
% mtabDetails(ks,:) = [];
% mtabNames(ks,:) = [];
% clear s ks

save(NameOfFile)

