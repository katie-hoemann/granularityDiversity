clear;
clc;

%% set parameters
numberSDs = 1.25; % number of SDs to use in thresholding out prompts

%% load and configure data file
rawData = importdata(['dataSet3_LIWC2015 results_150_10.csv']);
splitText1 = split(rawData.textdata(2:end,1),'_'); % split filename variable by _ to separate PPID and day number
splitText2 = split(splitText1(:,3),'.'); % split filename variable by . to separate out prompt ID
PPID = str2double(splitText1(:,1)); % convert to number
day = str2double(splitText1(:,2)); % convert to number
prompt = str2double(splitText2(:,1)); % convert to number
allData = [PPID day prompt rawData.data(:,2:end)]; % concatenate numeric data in single matrix
subjectIDlist = unique(allData(:,1)); % create subject ID list

%% load granularity and other derived variable(s)
load(['dataSet3_granularity_affect_measures.mat']); % load matrix of variables per subject (measures)
measures(:,2:4) = measures(:,2:4)*-1; % invert zICC values to get granularity (zInv)

%% grab data for each subject and run through calculations
themeMeans = mean(allData(:,4:end));
themeSDs = std(allData(:,4:end));
numThemes = size(allData(:,4:end),2);
for i_theme = 1:numThemes
    allData(:,3+i_theme) = allData(:,3+i_theme)>(themeMeans(i_theme)+numberSDs*themeSDs(i_theme));
    propDataTheme(i_theme,:) = sum(allData(:,3+i_theme))/size(allData,1);
end
mean(propDataTheme) % display average proportion of texts retained per theme
for i_subject = 1:length(subjectIDlist)
    subjectData = [];
    subjectID = subjectIDlist(i_subject);
    index = find(allData(:,1)==subjectID);
    subjectData = allData(index,4:end);
    % calculate number of prompts
    numPrompts(i_subject,1) = size(subjectData,1);
    % calculate averge number of themes per prompt, number prompts with no theme
    meanThemes(i_subject,1) = sum(subjectData>0,'all')/numPrompts(i_subject,1);
    noThemes(i_subject,1) = sum(sum(subjectData,2)==0)/numPrompts(i_subject,1);
    % calculate thematic diversity using Shannon formula (Quoidbach et al, 2014)
    totalCount(i_subject,1) = sum(sum(subjectData>0));
    for i_theme = 1:numThemes
        countTheme = sum(subjectData(:,i_theme)>0);
        propTheme(i_theme) = (countTheme./totalCount(i_subject,1))*log(countTheme./totalCount(i_subject,1));
    end
    propTheme = propTheme(isfinite(propTheme)); % remove values that are not finite (i.e., themes never expressed)
    diversityTheme(i_subject,1) = -1*sum(propTheme);
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

% drop subjects excluded from granularity analysis
subjectIDlist1 = array2table(subjectIDlist,'VariableNames',{'PPID1'});
subjectIDlist2 = array2table(measures(:,1),'VariableNames',{'PPID2'});
subjectIDlist_join = outerjoin(subjectIDlist1,subjectIDlist2,'LeftKeys','PPID1','RightKeys','PPID2');
subjectsToDrop = find(isnan(subjectIDlist_join.PPID2));

diversityTheme(subjectsToDrop) = [];
giniCoefTheme(subjectsToDrop) = [];
meanThemes(subjectsToDrop) = [];
noThemes(subjectsToDrop) = [];
numPrompts(subjectsToDrop) = [];
subjectIDlist(subjectsToDrop) = [];

% drop subjects excluded from text analysis
measures(find(isnan(subjectIDlist_join.PPID1)),:) = []; 

%% run multiple linear regression accounting for mean affect, number of prompts
% MODEL 1: mean negative affect, number of prompts, negative granularity
y1_S = zscore(giniCoefTheme); % z-score variables to get standardized beta as output
x1_S = [zscore(measures(:,8)) zscore(numPrompts) zscore(measures(:,2))]; 
mdl1_S = fitlm(x1_S,y1_S); % run the model once to get diagnostics
outliers1 = find(mdl1_S.Diagnostics.CooksDistance>3*mean(mdl1_S.Diagnostics.CooksDistance)); % use 3*M(Cook's distance) to identify outliers
mdl1_S = fitlm(x1_S,y1_S,'exclude',outliers1); % re-run model with outliers removed

y1_S_NO = y1_S; % create copies of variables with outliers removed
y1_S_NO(outliers1) = [];
negGran_S_NO = zscore(measures(:,2));
negGran_S_NO(outliers1,:) = [];
x1_S_NO = x1_S;
x1_S_NO(outliers1,:) = [];
x1_S_NO_padded = [ones(size(x1_S_NO,1),1) x1_S_NO]; 
[~,mdl1_S_CI,~,~,~] = regress(y1_S_NO,x1_S_NO_padded); % get 95% CI for beta coefficients

x1_S_R = [zscore(measures(:,8)) zscore(numPrompts)]; % get residuals for plotting
mdl1_S_R = fitlm(x1_S_R,y1_S);
y1_S_R = table2array(mdl1_S_R.Residuals(:,1)); 
y1_S_R_NO = y1_S_R;
y1_S_R_NO(outliers1,:) = [];

% MODEL 2: mean positive affect, number of prompts, positive granularity
y2_S = zscore(giniCoefTheme); % z-score variables to get standardized beta as output
x2_S = [zscore(measures(:,7)) zscore(numPrompts) zscore(measures(:,3))]; 
mdl2_S = fitlm(x2_S,y2_S); % run the model once to get diagnostics
outliers2 = find(mdl2_S.Diagnostics.CooksDistance>3*mean(mdl2_S.Diagnostics.CooksDistance)); % use 3*M(Cook's distance) to identify outliers
mdl2_S = fitlm(x2_S,y2_S,'exclude',outliers2); % re-run model with outliers removed

mdl2_X_S = [ones(size(measures,1),1) x2_S];
[~,mdl2_S_CI,~,~,~] = regress(y1_S,mdl2_X_S); % get 95% CI for beta coefficients
predictors_S_R2 = [zscore(measures(:,7)) zscore(numPrompts)]; 
mdl2_S_R = fitlm(predictors_S_R2,y1_S);
y_S_R2 = table2array(mdl2_S_R.Residuals(:,1)); % get residuals for plotting

y2_S_NO = y2_S; % create copies of variables with outliers removed
y2_S_NO(outliers2) = [];
posGran_S_NO = zscore(measures(:,3));
posGran_S_NO(outliers2,:) = [];
x2_S_NO = x2_S;
x2_S_NO(outliers2,:) = [];
x2_S_NO_padded = [ones(size(x2_S_NO,1),1) x2_S_NO]; 
[~,mdl2_S_CI,~,~,~] = regress(y2_S_NO,x2_S_NO_padded); % get 95% CI for beta coefficients

x2_S_R = [zscore(measures(:,7)) zscore(numPrompts)]; % get residuals for plotting
mdl2_S_R = fitlm(x2_S_R,y2_S);
y2_S_R = table2array(mdl2_S_R.Residuals(:,1)); 
y2_S_R_NO = y2_S_R;
y2_S_R_NO(outliers2,:) = [];

%% summary table
summary = [table2array(mdl1_S.Coefficients(end,1)) (table2array(mdl1_S.Coefficients(end,4))/2)...
    mdl1_S_CI(end,1) mdl1_S_CI(end,2) table2array(mdl1_S.Coefficients(end,3)) mdl1_S.DFE numel(outliers1);...
    table2array(mdl2_S.Coefficients(end,1)) (table2array(mdl2_S.Coefficients(end,4))/2)...
    mdl2_S_CI(end,1) mdl2_S_CI(end,2) table2array(mdl2_S.Coefficients(end,3)) mdl2_S.DFE numel(outliers2)];
summary_Table = array2table(summary,'VariableNames',{'B','p','CI_l','CI_u','t','df','outliers'},...
    'RowNames',{'negative granularity','positive granularity'});

%% create scatter plots of emotional granularity against residuals
subplot(1,2,1);
scatter1 = scatter(negGran_S_NO,y1_S_R_NO,[],rgb('SteelBlue'),'filled'); 
set(gca,'fontsize',14)
xlim([-3 3]);
ylim([-3 3]);
axis square;
xlabel('negative emotional granularity');
ylabel('experiential diversity (residuals)');
h1 = lsline;
h1.Color = 'k';
hold on;
scatter2 = scatter(zscore(measures(outliers1,2)),y1_S_R(outliers1),[],rgb('SteelBlue'));
hold off;

subplot(1,2,2);
scatter3 = scatter(posGran_S_NO,y2_S_R_NO,[],rgb('FireBrick'),'filled'); 
set(gca,'fontsize',14)
xlim([-3 3]);
ylim([-3 3]);
axis square;
xlabel('positive emotional granularity');
ylabel('experiential diversity (residuals)');
h1 = lsline;
h1.Color = 'k';
hold on;
scatter2 = scatter(zscore(measures(outliers2,3)),y2_S_R(outliers2),[],rgb('FireBrick'));
hold off;