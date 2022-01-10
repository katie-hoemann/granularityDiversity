clear;
clc;

%% specify dataset and parameters
dataSet = 1; % set to 1, 2, or 3
print = 0;

%% load data file, along with word file that includes raw norms
if dataSet == 1
    dataFile = 'dataSet1_EOD_rating data.xlsx'; 
    wordFile = 'dataSet1_words.csv'; 
elseif dataSet == 2
    dataFile = 'dataSet2_ESM_rating data.xlsx';
    wordFile = 'dataSet2_words.csv';
elseif dataSet == 3
    dataFile = 'dataSet3_ESM_rating data.xlsx';
    wordFile = 'dataSet3_words.csv';
end
rawData = importdata(dataFile);
allData = rawData.data;
subjectIDlist = unique(rawData.data(:,1)); % grab subject IDs from first column of data file
words = readtable(wordFile); 
wordList = rawData.colheaders(2:end)';  % grab sampled words from top row of data file

%% set parameters
if dataSet == 2
    startRatingsat1 = 1;
    maximumRating = 5; % set to highest valid scale value
elseif dataSet == 3
    startRatingsat1 = 0;
    maximumRating = 100;
else
    startRatingsat1 = 0;
    maximumRating = 6;
end

%% set valence and arousal categories for sampled words
for i_word = 1:height(words) % define valence categories
if words.Valence(i_word) > 5 % derived based on the database mean for Warriner et al (2013)
    valCat(i_word) = {'Positive'};
    positive(i_word) = 1;
    valence(i_word) = 1;
else
    valCat(i_word) = {'Negative'};
    positive(i_word) = 0;
    valence(i_word) = 2;
end
end 

for i_word = 1:height(words) % define arousal categories
if words.Arousal(i_word) > 4.6 % derived based on the sample mean for 88 PANAS-X terms in Warriner et al (2013)
    aroCat(i_word) = {'High'};
    high(i_word) = 1;
else
    aroCat(i_word) = {'Low'};
    high(i_word) = 0;
end
end 

words = [words valCat' aroCat']; % append table with category assignments
words.Properties.VariableNames(5:end) = {'ValCat' 'AroCat'}; % label new variables
labels = [positive' high']; % create matrix for logical indexing in ICC commands

%% grab data for each subject and run through calculations
for i_subject = 1:length(subjectIDlist)
    subjectData = [];
    subjectID = subjectIDlist(i_subject);
    index1 = find(allData(:,1)==subjectID);
    subjectData = allData(index1,2:end);
    % remove invalid values and missing data
    subjectData(subjectData > maximumRating) = NaN;
    missingData = isnan(subjectData); 
    missingData2 = any(missingData,2); 
    subjectData = subjectData(~missingData2,:); 
    % if necessary, rescale data to start at 0
    if startRatingsat1 == 1
        subjectData = subjectData-1;
    end
    % calculate number of events
    numEvents(i_subject) = size(subjectData,1);
    % calculate granularity (ICCs and zICCs)
    rawICC(i_subject,1) = ICC(subjectData(:,(labels(:,1)==0)),'A-k'); % negative valence ICC
    rawICC(i_subject,2) = ICC(subjectData(:,(labels(:,1)==1)),'A-k'); % positive valence ICC
    rawICC(i_subject,3) = (rawICC(i_subject,1)+rawICC(i_subject,2))/2; % valence mean ICC
    rawICC(rawICC<0) = 0;
    rtoZ(i_subject,1) = 0.5*log((1+rawICC(i_subject,1))/(1-rawICC(i_subject,1)));
    rtoZ(i_subject,2) = 0.5*log((1+rawICC(i_subject,2))/(1-rawICC(i_subject,2)));
    rtoZ(i_subject,3) = 0.5*log((1+rawICC(i_subject,3))/(1-rawICC(i_subject,3)));
    % calculate emodiversity using Shannon formula (Quoidbach et al, 2014)
    totalCount = sum(sum(subjectData>0));
    for i_emotion = 1:length(wordList)
        countEmotion = sum(subjectData(:,i_emotion)>0);
        propEmotion(i_emotion) = (countEmotion./totalCount)*log(countEmotion./totalCount);
    end
    propEmotion = propEmotion(isfinite(propEmotion)); % remove values that are not finite (i.e., emotions never experienced)
    emodiversity(i_subject) = -1*sum(propEmotion);
    % calculate emodiversity using Gini coefficient (Benson et al, 2018)
    numEmotions = length(wordList);
    for i_emotion = 1:numEmotions
        countEmotionRanked(i_emotion) = sum(subjectData(:,i_emotion)>0);
    end
    countEmotionRanked = sort(countEmotionRanked);
    index = 1:1:numEmotions;
    for i_emotion = 1:numEmotions
        weightedCount(i_emotion) = countEmotionRanked(i_emotion)*index(i_emotion);
    end
    giniCoefEmo(i_subject) = 1-(((2*sum(weightedCount))/(numEmotions*sum(countEmotionRanked)))-((numEmotions+1)/numEmotions));
    % calculate mean and SD for positive and negative valence
    mPositive(i_subject) = mean(mean(subjectData(:,(valence==1))));
    mNegative(i_subject) = mean(mean(subjectData(:,(valence==2))));
    sdPositive(i_subject) = mean(std(subjectData(:,(valence==1))));
    sdNegative(i_subject) = mean(std(subjectData(:,(valence==2)))); 
end

%% create summary table and save data
measures = horzcat(subjectIDlist, rtoZ, emodiversity', giniCoefEmo', mPositive', mNegative', sdPositive', sdNegative', numEvents');
save(['dataSet' num2str(dataSet) '_granularity_affect_measures.mat'],'measures');
variableNames = {'PPID','zICC_N','zICC_P','zICC_M','emodiv','giniEmo','mPos','mNeg','sdPos','sdNeg','numEvents'};
measures_Table = array2table(measures,'VariableNames',variableNames);