%get the MAVEN data and the sample information squared away
%working on the Chisholm samples run on TSQ, pos and neg different
%allow poor quality scores for the new ''methyl-*tanoic' acids
%KL 10/20/2016; KL 11/14/2017 add in the ability to merge based on names
%KL 1/20/2018 - around line 245 is code that works here (MATLAB 2013b, but
%will not work in MATLAB 2017a). 
%KL 12/11/2018 - working on year 3 of BIOS-SCOPE data
%KL 2/21/2020 working on Pro8, define SRM files here 
%KL 3/16/2021 updating to use Melissa's new names
%KL 4/10/2021 - updated to use three runs on TSQ
clear

%%%set the file names up front:
%MAT file for MATLAB; CHANGE NAME
NameOfFile = 'TSQ_KHU8.2021.06.28.mat';

%sequence file from the TSQ; CHANGE PATH
wDir = '/Users/kathrynhalloran/Documents/MATLAB/Uptake experiments/KHU8';
fName = 'mtab_KHU8_2021_0618_edit.xls';

sampleInfoFile = [wDir filesep fName];
%Since there is no comma at the end of the headerline, MATLAB will crash
%with the import...read in the first line, make sure it ends in a comma, 
%close the file and then open it again...

%%SRM files; CHANGE PATH
sDir = '/Users/kathrynhalloran/Documents/MATLAB/Uptake experiments/KHU8';
%SRM_pos = [sDir filesep 'SRMList_Pos_022521_UptakeIncubations_forElMAVEN_v1.csv'];
%SRM_neg = [sDir filesep 'SRMList_Neg_022521_UptakeIncubations_forElMAVEN_v4.csv'];
SRM_mix8 = [sDir filesep 'mtab_KHU8_2021_0617_SRMList_edit.csv'];
clear sDir

%since I always seem to have part1 and part2 files from MAVEN (bc split
%with manual vs. automatic integration, make MATLAB combine the files
%note that I now have 8 files (6 in positive ion mode and 2 in negative)
%readFile1_pos = 'TSQ_EE2_pos_part1.2021.04.09.csv';
%readFile2_pos = 'TSQ_EE2_pos_part2.2021.04.09.csv';

%readFile1_neg = 'TSQ_EE2_neg_part1.2021.04.09.csv';
%readFile2_neg = 'TSQ_EE2_neg_part2.2021.04.09.csv';
%readFile3_neg = 'TSQ_EE2_neg_part3.2021.04.09.csv';
%readFile4_neg = 'TSQ_EE2_neg_part4.2021.04.09.csv';

readFile1_mix8 = 'KHU8_peaklist_edit.csv';

%setup some variables for the TSQ processing:
setQuality = 0.1;
warning('off', 'MATLAB:codetools:ModifiedVarnames');

%%first deal with the positive ion mode files
%%first deal with the positive ion mode files
% tempFile = 'temp.csv';
% fid1 = fopen(readFile1_pos);
% tline = fgetl(fid1);
% if ~strcmp(tline(end),',')
%     tline = strcat(tline,','); %put the comma at the end of the headerline
% end
% %open up the temp file, and set it to be write-able
% fidOut = fopen(tempFile,'w');
% fprintf(fidOut, '%s\n',tline);
% 
% tline = fgetl(fid1); %get the next line of the file
% while ischar(tline)
%     fprintf(fidOut, '%s\n', tline);
%     tline = fgetl(fid1);    
% end
% fclose(fid1);
% 
% %now get the second file...since the headerline will be the same can skip it
% fid2 = fopen(readFile2_pos);
% tline = fgetl(fid2);
% tline = fgetl(fid2); %get the next line of the file
% while ischar(tline)
%     fprintf(fidOut, '%s\n', tline);
%     tline = fgetl(fid2);    
% end
% fclose(fid2);
% 
% % %for rest of files files...
% % fid3 = fopen(readFile3_pos);
% % tline = fgetl(fid3);
% % tline = fgetl(fid3); %get the next line of the file
% % while ischar(tline)
% %     fprintf(fidOut, '%s\n', tline);
% %     tline = fgetl(fid3);    
% % end
% % fclose(fid3);   
% % 
% % %for rest of files files...
% % fid4 = fopen(readFile4_pos);
% % tline = fgetl(fid4);
% % tline = fgetl(fid4); %get the next line of the file
% % while ischar(tline)
% % fprintf(fidOut, '%s\n', tline);
% % tline = fgetl(fid4);    
% % end
% % fclose(fid4);   
%        
% 
% fclose(fidOut);
% clear fid* tline fidOut ans readFile*_pos 
% 
% %%%now move on with the TSQ data processing
% [pos.sNames pos.kgd pos.curves] = considerMAVENv14_KHU8(tempFile,sampleInfoFile,setQuality,1,'positive',SRM_pos);
% %there are zeros and NaNs in the kgd.goodData...
% i = isnan(pos.kgd.goodData);
% pos.kgd.goodData(i) = 0;
% clear i
% 
% 
% %%now do the negative ion mode files
% %%now do the negative ion mode files
% tempFile = 'temp.csv';
% fid1 = fopen(readFile1_neg);
% tline = fgetl(fid1);
% if ~strcmp(tline(end),',')
%     tline = strcat(tline,','); %put the comma at the end of the headerline
% end
% %open up the temp file, and set it to be write-able
% fidOut = fopen(tempFile,'w');
% fprintf(fidOut, '%s\n',tline);
% 
% tline = fgetl(fid1); %get the next line of the file
% while ischar(tline)
%     fprintf(fidOut, '%s\n', tline);
%     tline = fgetl(fid1);    
% end
% fclose(fid1);
% 
% %now get the second file...since the headerline will be the same can skip it
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
%     %for three files...
%     fid3 = fopen(readFile3_neg);
%     tline = fgetl(fid3);
%     tline = fgetl(fid3); %get the next line of the file
%     while ischar(tline)
%         fprintf(fidOut, '%s\n', tline);
%         tline = fgetl(fid3);    
%     end
%     fclose(fid3);    
% end
% 
% if 1
%     %for four files...
%     fid4 = fopen(readFile4_neg);
%     tline = fgetl(fid4);
%     tline = fgetl(fid4); %get the next line of the file
%     while ischar(tline)
%         fprintf(fidOut, '%s\n', tline);
%         tline = fgetl(fid4);    
%     end
%     fclose(fid4);    
% end
% 
% fclose(fidOut);
% clear fid* tline fidOut ans readFile*_neg 
% 
% %%%now move on with the TSQ data processing
% [neg.sNames neg.kgd neg.curves] = considerMAVENv14_KHU8(tempFile,sampleInfoFile,setQuality,1,'negative',SRM_neg);
% %housecleaning
% i = isnan(neg.kgd.goodData);
% neg.kgd.goodData(i) = 0;
% clear i
% clear wDir fName tempFile 


%%now do the Mix 8 files
%%now do the Mix 8 files
tempFile = 'temp.csv';
fid1 = fopen(readFile1_mix8);
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
fclose(fidOut);

clear fid1 fid2 fid3 tline fidOut ans readFile1_neg readFile2_neg readFile3_neg

%%%now move on with the TSQ data processing
[both.sNames both.kgd both.curves] = considerMAVENv14_KHU8(tempFile,sampleInfoFile,setQuality,1,'both',SRM_mix8);
%housecleaning
i = isnan(both.kgd.goodData);
both.kgd.goodData(i) = 0;
clear i
clear wDir fName tempFile 
clear setQuality setErrorP
clear SRM_*
clear readFile*_*


%now, I have three structures with data (pos, neg, mix8) - merge all three

%get all the metabolite names...should be a unique list even when I merge
%across positive and negative ion mode, but check that

%as expected, there are issues. For all the mix 8 standards, I need a
% %marker so I know this is mix8
for a = 1:size(both.kgd.name,1)
    both.kgd.name(a,1) = strcat(both.kgd.name(a),{' '},'mix8');
end
clear a

mtabNames = sort(cat(1,both.kgd.name));
if length(unique(mtabNames)) ~= length(mtabNames)
    error('Something is wrong - duplicate names in the list of metabolites')
end

%for the BIOS-SCOPE temporal data, I will have duplicate sets of names with
%either _pos or _neg appended; they are in the first worksheet
tInfo = readtable(sampleInfoFile);
clear sampleInfoFile

%first, go through and iterate through the pooled samples and monster
%samples to provide numbers for these (otherwise will have duplicate
%names)
%% only need to do this if the pooled samples are not numbered
%%
if 0
    s = strcmp(tInfo.SampleName,'EE_Pool pos');
    ks = find(s==1);
    for a = 1:length(ks);
        t = tInfo.SampleName(ks(a));
        tInfo.SampleName(ks(a)) = strcat('p',num2str(a),t);
        clear t
    end
    clear a ks a

    s = strcmp(tInfo.SampleName,'EE_pool neg');
    ks = find(s==1);
    for a = 1:length(ks);
        t = tInfo.SampleName(ks(a));
        tInfo.SampleName(ks(a)) = strcat('p',num2str(a),t);
        clear t
    end
    clear a ks a   

    s = strcmp(tInfo.SampleName,'EE_pool mix8');
    ks = find(s==1);
    for a = 1:length(ks);
        t = tInfo.SampleName(ks(a));
        tInfo.SampleName(ks(a)) = strcat('p',num2str(a),t);
        clear t
    end
    clear a ks a   
end


%before I dive into the unknowns, remove anything that has goodData = 0
%at this point, also only want Unknowns
s = strcmp(tInfo.SampleType,'Unknown');
k = find(tInfo.goodData==0 | s ~= 1);
tInfo(k,:) = [];
clear k s

%now find the Unknown...should have the same number for positive and
%negative ion mode bc have pruned out the different QC samples already
s = strcmp(tInfo.SampleType,'Unknown');
sp = strcmp(tInfo.ionMode,'positive');
ksp = (find(s==1 & sp==1));
sn = strcmp(tInfo.ionMode,'negative');
ksn = (find(s==1 & sn==1));
sb = strcmp(tInfo.ionMode,'both');
ksb = (find(s==1 & sn==1));


if ~isequal(length(ksp),length(ksn),length(ksb))
    error('Something wrong, these should be the same length')
end
clear a sp sn sb ksp ksn ksb


sInfo = table;
sInfo.cName = unique(tInfo.shortName);
%the first row of this will be empty, delete that
if isequal(sInfo.cName(1),{''});
    sInfo(1,:) = [];
end


%now make an empty matrix for the data...will be all numbers so no need for
%special format
mtabData = zeros(size(mtabNames,1),size(sInfo,1));
%need to track some additional details:
mtabDetails = table();

%get the index for rows for positive AND negative mtabs:
%[c idx_posNew idx_posOld] = intersect(mtabNames,pos.kgd.name);
%[c idx_negNew idx_negOld] = intersect(mtabNames,neg.kgd.name);
[c idx_bothNew idx_bothOld] = intersect(mtabNames,both.kgd.name);

%mtabDetails.mode(idx_posNew,1) = {'pos'};
%mtabDetails.mode(idx_negNew,1) = {'neg'};
mtabDetails.mode(idx_bothNew,1) = {'both'};

%sInfo.runOrder_pos(:,1) = 0;
%sInfo.runOrder_neg(:,1) = 0;
sInfo.runOrder_both(:,1) = 0;

%sInfo.FileName_pos(:,1) = {''};
%sInfo.FileName_neg(:,1) = {''};
sInfo.FileName_both(:,1) = {''};


% for a = 1:size(sInfo,1);
%     s = strcmp(sInfo.cName(a),tInfo.shortName);
%     ks = find(s==1);
%     if length(ks) ~= 3
%         sInfo.cName(a)
%         error('Something is wrong, should be three of each')
%     end
% %     
% %     %some variant of this:
% %     for aa = 1:3
% %         %propagate sInfo with the basic information, only do once
% %         %no need to completely fill this in bc not working on data analysis
% %         if 1
% %             if aa == 1
% %             
% %             sInfo.type(a) = tInfo.type(ks(aa));
% %             sInfo.replicate(a) = tInfo.replicate(ks(aa));
% % %             sInfo.otherName(a) = tInfo.otherName(ks(aa));
% % %             sInfo.volSent_ml(a) = tInfo.volSent_ml(ks(aa));
% % %             sInfo.cellsML(a) = tInfo.cellsML(ks(aa));
% % %             sInfo.cellsFiltered(a) = tInfo.cellsFiltered(ks(aa));
% %             end
% %         end
% %         
% %         im = tInfo.ionMode{ks(aa)};
% %         if isequal(im,'positive')
% %             tName = tInfo.FileName(ks(aa));
% %             sInfo.FileName_pos(a,1) = tName;
% % 
% %             [c ia tIdx] =intersect(tName,pos.sNames);
% %             mtabData(idx_posNew,a) = pos.kgd.goodData(idx_posOld,tIdx);
% %             clear c ia tIdx tName
% %             
% %         elseif isequal(im,'negative')
% %             tName = tInfo.FileName(ks(aa));
% %             sInfo.FileName_neg(a,1) = tName;
% % 
% %             [c ia tIdx] =intersect(tName,neg.sNames);
% %             mtabData(idx_negNew,a) = neg.kgd.goodData(idx_negOld,tIdx);
% %             clear c ia tIdx tName
% % 
% %         elseif isequal(im,'both')
% %             tName = tInfo.FileName(ks(aa));
% %             sInfo.FileName_both(a,1) = tName;
% % 
% %             [c ia tIdx] =intersect(tName,both.sNames);
% %             mtabData(idx_bothNew,a) = both.kgd.goodData(idx_bothOld,tIdx);
% %             clear c ia tIdx tName
% %             
% %         else 
% %             error('Something wrong')
% %         end
% %         clear im
% %     end
% %     clear aa s ks        
% end
clear a

clear idx_* tInfo

for a = 1: size(sInfo,1)
    %do positive ion mode first
%     gc = sInfo{a,'FileName_pos'}{:}; %added {:} to deal with table output
%     t = regexp(gc,'_');
%     if ~isempty(t)
%         sInfo.runOrder_pos(a,1) = str2num(gc(t(end)+1:end));
%     else
%         sInfo.runOrder_pos(a,1) = NaN;
%     end
%     clear gc t
%     
%     %then negative ion mode 
%     gc = sInfo{a,'FileName_neg'}{:}; %added {:} to deal with table output
%     t = regexp(gc,'_');
%     if ~isempty(t)
%         sInfo.runOrder_neg(a,1) = str2num(gc(t(end)+1:end));
%     else
%         sInfo.runOrder_neg(a,1) = NaN;
%     end
%     clear gc t
    
    %then Mix8 (both: positive and negative ion mode) 
    gc = sInfo{a,'FileName_both'}{:}; %added {:} to deal with table output
    t = regexp(gc,'_');
    if ~isempty(t)
        sInfo.runOrder_both(a,1) = str2num(gc(t(end)+1:end));
    else
        sInfo.runOrder_both(a,1) = NaN;
    end
    clear gc t
end
clear a
 
if 0
    % %making decisions about pos/neg mtabs (3/15/2021 not updated yet)
    toDelete = {'NAD neg','adenine neg','adenosine 5''-monophosphate neg'...
        'desthiobiotin neg','guanosine neg',...
        'hemin I','n-acetyl muramic acid neg','pantothenic acid pos',...
        's-(5''-adenosyl)-L-homocysteine neg', 'xanthine neg',...
        'xanthosine pos','taurocholic acid d5',...
        'glutamic acid d3 pos','cytidine neg',...
        'sn-glycerol 3-phosphate pos ','chitobiose pos',...
        'biotin neg','folic acid neg','meso26diaminopimelicacid',...
        '2deoxyinosine neg','aspartic acid pos'};
    [c ia ib] = intersect(toDelete,mtabNames);
    
    if ~isequal(length(c),length(toDelete))
        missing = setdiff(toDelete,mtabNames)
        error('Something is wrong, these should be the same length')
    end

    mtabNames(ib,:)=[];
    mtabDetails(ib,:)=[];
    mtabData(ib,:)=[];
    clear c ia ib
    clear toDelete
    
end
% 
% %for glutathione oxidized, need to set the filter samples to zeros
% fprintf('Glutathione oxidized in filters set to zero\n')
% sr = strcmp(mtabNames,'l-glutathione oxidized pos');
% ksr = find(sr==1);
% sc = strcmp(sInfo.sampleSource,'filter');
% ksc = find(sc==1);
% 
% mtabData(ksr,ksc) = 0;
% clear sr ksr sc ksc


save(NameOfFile)
