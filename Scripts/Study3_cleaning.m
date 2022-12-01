clear all;
clc;

%% import data
% event descriptions/texts
[~, ~, rawTextData] = xlsread('dataSet3_ESM_event descriptions.csv'); % import the data in its original format (i.e., both string and numerical)
numericData = rawTextData(2:end,[1:4 6]);
textData = rawTextData(2:end,5);
subjectIDlistText = unique(cell2mat(numericData(:,2)));

% emotion intensity ratings
rawRatingData = importdata('dataSet3_ESM_rating data.xlsx');
subjectIDlistRating = unique(rawRatingData.data(:,1));

%% find & remove excluded participants
% identify participants to drop
subjectIDlist1 = array2table(subjectIDlistText,'VariableNames',{'PPID1'});
subjectIDlist2 = array2table(subjectIDlistRating,'VariableNames',{'PPID2'});
subjectIDlist_join = outerjoin(subjectIDlist1,subjectIDlist2,'LeftKeys','PPID1','RightKeys','PPID2');
subjectIndicesToDrop = find(isnan(subjectIDlist_join.PPID2));
subjectsToDrop = subjectIDlistText(subjectIndicesToDrop);

% remove participant data from text file
for i_subject = 1:length(subjectsToDrop)
    index = cell2mat(numericData(:,2)) == subjectsToDrop(i_subject); 
    numericData(index,:) = [];
    textData(index) = [];
end
newSubjectIDlistText = unique(cell2mat(numericData(:,2)));
lengthText = cell2mat(numericData(:,5));

%% get descriptive stats & plot data
minLength = min(lengthText);
maxLength = max(lengthText);
mdnLength = median(lengthText);
mLength = mean(lengthText);
sdLength = std(lengthText);
histogram(lengthText);
% participant-level stats
for i_subject = 1:length(newSubjectIDlistText)
    subjectID = newSubjectIDlistText(i_subject);
    index1 = find(cell2mat(numericData(:,2))==subjectID);
    numEntries(i_subject,:) = length(cell2mat(numericData(index1,4)));
end

% try removing texts < 30 words long
textData_30 = textData;
numericData_30 = numericData;
textData_30(lengthText<30) = [];
numericData_30(lengthText<30,:) = [];
percentRemoved_30 = (length(textData)-length(textData_30))/length(textData);
mdnLength_30 = median(lengthText(lengthText>=30));
mLength_30 = mean(lengthText(lengthText>=30));
sdLength_30 = std(lengthText(lengthText>=30));
% participant-level stats
for i_subject = 1:length(newSubjectIDlistText)
    subjectID = newSubjectIDlistText(i_subject);
    index1 = find(cell2mat(numericData_30(:,2))==subjectID);
    numEntries_30(i_subject,:) = length(cell2mat(numericData_30(index1,4)));
end
loss_30 = (numEntries-numEntries_30)./numEntries*100;
histogram(loss_30,'BinWidth',5);

% try removing texts < 25 words long
textData_25 = textData;
numericData_25 = numericData;
textData_25(lengthText<25) = [];
numericData_25(lengthText<25,:) = [];
percentRemoved_25 = (length(textData)-length(textData_25))/length(textData);
mdnLength_25 = median(lengthText(lengthText>=25));
mLength_25 = mean(lengthText(lengthText>=25));
sdLength_25 = std(lengthText(lengthText>=25));
% participant-level stats
for i_subject = 1:length(newSubjectIDlistText)
    subjectID = newSubjectIDlistText(i_subject);
    index1 = find(cell2mat(numericData_25(:,2))==subjectID);
    numEntries_25(i_subject,:) = length(cell2mat(numericData_25(index1,4)));
end
loss_25 = (numEntries-numEntries_25)./numEntries*100;
histogram(loss_25,'BinWidth',5);

% try removing texts < 20 words long
textData_20 = textData;
numericData_20 = numericData;
textData_20(lengthText<20) = [];
numericData_20(lengthText<20,:) = [];
percentRemoved_20 = (length(textData)-length(textData_20))/length(textData);
mdnLength_20 = median(lengthText(lengthText>=20));
mLength_20 = mean(lengthText(lengthText>=20));
sdLength_20 = std(lengthText(lengthText>=20));
% participant-level stats
for i_subject = 1:length(newSubjectIDlistText)
    subjectID = newSubjectIDlistText(i_subject);
    index1 = find(cell2mat(numericData_20(:,2))==subjectID);
    numEntries_20(i_subject,:) = length(cell2mat(numericData_20(index1,4)));
end
loss_20 = (numEntries-numEntries_20)./numEntries*100;
histogram(loss_20,'BinWidth',5);

%% export data
newData = [numericData(:,1:4) textData];
newData = [{'Filename','PPID','Date','Time','Text'}; newData]; % optional headers
xlswrite('dataSet3_ESM_event descriptions_final N.csv',newData);

newData_30 = [numericData_30(:,1:4) textData_30];
newData_30 = [{'Filename','PPID','Date','Time','Text'}; newData_30]; % optional headers
xlswrite('dataSet3_ESM_event descriptions_final N_min 30.csv',newData_30);

newData_25 = [numericData_25(:,1:4) textData_25];
newData_25 = [{'Filename','PPID','Date','Time','Text'}; newData_25]; % optional headers
xlswrite('dataSet3_ESM_event descriptions_final N_min 25.csv',newData_25);

newData_20 = [numericData_20(:,1:4) textData_20];
newData_20 = [{'Filename','PPID','Date','Time','Text'}; newData_20]; % optional headers
xlswrite('dataSet3_ESM_event descriptions_final N_min 20.csv',newData_20);