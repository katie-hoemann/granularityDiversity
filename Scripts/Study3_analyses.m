clear;
clc;

%% set parameters
threshold = 1; % set to 1 to threshold prompts based on SDs
numberSDs = 1.25; % number of SDs to use in thresholding out prompts
typeICC = 'C-k'; % set to A-k for agreement; C-k for consistency
generateFigures = 1; % set to 1 to generate summary figures

%% load and configure data
% LIWC data (experiential diversity)
rawLIWCData = importdata(['dataSet3_min25_LIWC2015 results_150_10_min25.csv']);
splitText1 = split(rawLIWCData.textdata(2:end,1),'_'); % split filename variable by _ to separate PPID and day number
splitText2 = split(splitText1(:,3),'.'); % split filename variable by . to separate out prompt ID
PPID = str2double(splitText1(:,1)); % convert to number
day = str2double(splitText1(:,2)); % convert to number
prompt = str2double(splitText2(:,1)); % convert to number
allLIWCData = [PPID day prompt rawLIWCData.data(:,2:end)]; % concatenate numeric data in single matrix
subjectIDlistLIWC = unique(allLIWCData(:,1)); % create subject ID list

% rating data (emotional granularity)
ratingDataFile = 'dataSet3_ESM_rating data.xlsx'; 
wordFile = 'dataSet3_words.csv'; 
rawRatingData = importdata(ratingDataFile);
allRatingData = rawRatingData.data;
subjectIDlistRating = unique(rawRatingData.data(:,1)); % grab subject IDs from first column of data file
words = readtable(wordFile); 
wordList = rawRatingData.colheaders(2:end)';  % grab sampled words from top row of data file

%% grab data for each subject and run through calculations
% experiential diversity
% first, threshold data and get the average proportion of texts scored above threshold
themeMeans = mean(allLIWCData(:,4:end));
themeSDs = std(allLIWCData(:,4:end));
numThemes = size(allLIWCData(:,4:end),2);
for i_theme = 1:numThemes
    if threshold == 1
        allLIWCData(:,3+i_theme) = allLIWCData(:,3+i_theme)>(themeMeans(i_theme)+numberSDs*themeSDs(i_theme));
    end
    propDataTheme(i_theme,:) = sum(allLIWCData(:,3+i_theme))/size(allLIWCData,1);
end
mean(propDataTheme) % display

for i_subject = 1:length(subjectIDlistLIWC)
    subjectData = [];
    subjectID = subjectIDlistLIWC(i_subject);
    index = find(allLIWCData(:,1)==subjectID);
    subjectData = allLIWCData(index,4:end);
    % calculate descriptive statistics
    numPrompts(i_subject,1) = size(subjectData,1);
    totalCount(i_subject,1) = sum(sum(subjectData>0));
    % calculate thematic diversity using Gini coefficient (Benson et al, 2018)
    for i_theme = 1:numThemes
        countThemeRanked(i_theme) = sum(subjectData(:,i_theme)>0);
    end
    countThemeRanked = sort(countThemeRanked);
    index = 1:1:numThemes;
    for i_theme = 1:numThemes
        weightedCount(i_theme) = countThemeRanked(i_theme)*index(i_theme);
    end
    giniCoefTheme(i_subject,1) = 1-(((2*sum(weightedCount))/(numThemes*sum(countThemeRanked)))-((numThemes+1)/numThemes));
    if isnan(giniCoefTheme(i_subject,1)) == 1
        giniCoefTheme(i_subject,1) = 0;
    end
end

% emotional granularity
% first, set valence and arousal categories for sampled words
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

for i_subject = 1:length(subjectIDlistRating)
    subjectData = [];
    subjectID = subjectIDlistRating(i_subject);
    index1 = find(allRatingData(:,1)==subjectID);
    subjectData = allRatingData(index1,2:end);
    % remove invalid values and missing data
    subjectData(subjectData > 100) = NaN;
    missingData = isnan(subjectData); 
    missingData2 = any(missingData,2); 
    subjectData = subjectData(~missingData2,:); 
    % calculate number of events
    numEvents(i_subject) = size(subjectData,1);
    % calculate granularity (ICCs and zICCs)
    rawICC(i_subject,1) = ICC(subjectData(:,(labels(:,1)==0)),typeICC); % negative valence ICC
    rawICC(i_subject,2) = ICC(subjectData(:,(labels(:,1)==1)),typeICC); % positive valence ICC
    rawICC(rawICC<0) = 0;
    rtoZ(i_subject,1) = 0.5*log((1+rawICC(i_subject,1))/(1-rawICC(i_subject,1)));
    rtoZ(i_subject,2) = 0.5*log((1+rawICC(i_subject,2))/(1-rawICC(i_subject,2)));
    negGran(i_subject,1) = rtoZ(i_subject,1)*-1; % invert zICC values to get intuitively scaled negative granularity
    posGran(i_subject,1) = rtoZ(i_subject,2)*-1; % invert zICC values to get intuitively scaled positive granularity
    % calculate mean positive and negative valence
    mPositive(i_subject,1) = mean(mean(subjectData(:,(valence==1))));
    mNegative(i_subject,1) = mean(mean(subjectData(:,(valence==2))));
end

% drop subjects excluded from text analysis
subjectIDlist1 = array2table(subjectIDlistLIWC,'VariableNames',{'PPID1'});
subjectIDlist2 = array2table(subjectIDlistRating,'VariableNames',{'PPID2'});
subjectIDlist_join = outerjoin(subjectIDlist1,subjectIDlist2,'LeftKeys','PPID1','RightKeys','PPID2');
subjectsToDrop = find(isnan(subjectIDlist_join.PPID1));

negGran(subjectsToDrop) = [];
posGran(subjectsToDrop) = [];
mNegative(subjectsToDrop) = [];
mPositive(subjectsToDrop) = [];
subjectIDlistRating(subjectsToDrop) = [];

%% run multiple linear regressions
% MODEL 1: predicting negative granularity from mean negative affect, number of prompts, experiential diversity
y1_S = zscore(negGran); % z-score variables to get standardized beta as output
x1_S = [zscore(mNegative) zscore(numPrompts) zscore(giniCoefTheme)]; 
mdl1_S = fitlm(x1_S,y1_S); % run the model once to get diagnostics
outliers1 = find(mdl1_S.Diagnostics.CooksDistance>3*mean(mdl1_S.Diagnostics.CooksDistance)); % use 3*M(Cook's distance) to identify outliers
mdl1_S = fitlm(x1_S,y1_S,'exclude',outliers1); % re-run model with outliers removed

y1_S_NO = y1_S; % create copies of variables with outliers removed
y1_S_NO(outliers1) = [];
giniCoefTheme_negGran_S_NO = zscore(giniCoefTheme);
giniCoefTheme_negGran_S_NO(outliers1,:) = [];
x1_S_NO = x1_S;
x1_S_NO(outliers1,:) = [];
x1_S_NO_padded = [ones(size(x1_S_NO,1),1) x1_S_NO]; 
[~,mdl1_S_CI,~,~,~] = regress(y1_S_NO,x1_S_NO_padded); % get 95% CI for beta coefficients

x1_S_R = [zscore(mNegative) zscore(numPrompts)]; % get residuals for plotting
mdl1_S_R = fitlm(x1_S_R,y1_S);
y1_S_R = table2array(mdl1_S_R.Residuals(:,1)); 
y1_S_R_NO = y1_S_R;
y1_S_R_NO(outliers1,:) = [];

% MODEL 2: predicting positive granularity from mean positive affect, number of prompts, experiential diversity
y2_S = zscore(posGran); % z-score variables to get standardized beta as output
x2_S = [zscore(mPositive) zscore(numPrompts) zscore(giniCoefTheme)]; 
mdl2_S = fitlm(x2_S,y2_S); % run the model once to get diagnostics
outliers2 = find(mdl2_S.Diagnostics.CooksDistance>3*mean(mdl2_S.Diagnostics.CooksDistance)); % use 3*M(Cook's distance) to identify outliers
mdl2_S = fitlm(x2_S,y2_S,'exclude',outliers2); % re-run model with outliers removed

y2_S_NO = y2_S; % create copies of variables with outliers removed
y2_S_NO(outliers2) = [];
giniCoefTheme_posGran_S_NO = zscore(giniCoefTheme);
giniCoefTheme_posGran_S_NO(outliers2,:) = [];
x2_S_NO = x2_S;
x2_S_NO(outliers2,:) = [];
x2_S_NO_padded = [ones(size(x2_S_NO,1),1) x2_S_NO]; 
[~,mdl2_S_CI,~,~,~] = regress(y2_S_NO,x2_S_NO_padded); % get 95% CI for beta coefficients

x2_S_R = [zscore(mPositive) zscore(numPrompts)]; % get residuals for plotting
mdl2_S_R = fitlm(x2_S_R,y2_S);
y2_S_R = table2array(mdl2_S_R.Residuals(:,1)); 
y2_S_R_NO = y2_S_R;
y2_S_R_NO(outliers2,:) = [];

%% create summary tables
summary_granularity = [table2array(mdl1_S.Coefficients(end,1)) (table2array(mdl1_S.Coefficients(end,4))/2)...
    mdl1_S_CI(end,1) mdl1_S_CI(end,2) table2array(mdl1_S.Coefficients(end,3)) mdl1_S.DFE numel(outliers1);...
    table2array(mdl2_S.Coefficients(end,1)) (table2array(mdl2_S.Coefficients(end,4))/2)...
    mdl2_S_CI(end,1) mdl2_S_CI(end,2) table2array(mdl2_S.Coefficients(end,3)) mdl2_S.DFE numel(outliers2)];
summary_granularity_Table = array2table(summary_granularity,'VariableNames',{'B','p','CI_l','CI_u','t','df','outliers'},...
    'RowNames',{'negative granularity','positive granularity'});

summary_negGran_model = [table2array(mdl1_S.Coefficients(:,1)) mdl1_S_CI(:,1) mdl1_S_CI(:,2) table2array(mdl1_S.Coefficients(:,3))...
    repmat(mdl1_S.DFE,mdl1_S.NumCoefficients,1) table2array(mdl1_S.Coefficients(:,4))/2 repmat(numel(outliers1),mdl1_S.NumCoefficients,1)];
summary_negGran_model_Table = array2table(summary_negGran_model,'VariableNames',{'B','CI_l','CI_u','t','df','p','outliers'},...
    'RowNames',{'intercept','mean negative affect','number of texts','experiential diversity'});

summary_posGran_model = [table2array(mdl2_S.Coefficients(:,1)) mdl2_S_CI(:,1) mdl2_S_CI(:,2) table2array(mdl2_S.Coefficients(:,3))...
    repmat(mdl2_S.DFE,mdl2_S.NumCoefficients,1) table2array(mdl2_S.Coefficients(:,4))/2 repmat(numel(outliers2),mdl2_S.NumCoefficients,1)];
summary_posGran_model_Table = array2table(summary_posGran_model,'VariableNames',{'B','CI_l','CI_u','t','df','p','outliers'},...
    'RowNames',{'intercept','mean positive affect','number of texts','experiential diversity'});

%% create scatter plots of emotional granularity against residuals
if generateFigures == 1
    s = 100;

    figure;
    scatter1 = scatter(giniCoefTheme_negGran_S_NO,y1_S_R_NO,s,rgb('SteelBlue'),'filled'); 
    set(gca,'fontsize',20)
    xlim([-3 3]);
    ylim([-3 3]);
    axis square;
    xlabel('experiential diversity');
    ylabel('negative emotional granularity (residuals)');
    h1 = lsline;
    h1.Color = 'k';
    hold on;
    scatter2 = scatter(zscore(giniCoefTheme(outliers1)),y1_S_R(outliers1),s,rgb('SteelBlue'));
    hold off;

    figure;
    scatter3 = scatter(giniCoefTheme_posGran_S_NO,y2_S_R_NO,s,rgb('FireBrick'),'filled'); 
    set(gca,'fontsize',20)
    xlim([-3 3]);
    ylim([-3 3]);
    axis square;
    xlabel('experiential diversity');
    ylabel('positive emotional granularity (residuals)');
    h1 = lsline;
    h1.Color = 'k';
    hold on;
    scatter2 = scatter(zscore(giniCoefTheme(outliers2)),y2_S_R(outliers2),s,rgb('FireBrick'));
    hold off;
end