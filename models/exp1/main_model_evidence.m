clear variables
close all

load basic_cvll.mat

%% Calculate
subInd = [1:24];

fnames = fieldnames(cvll_basic);

n = 0;
for s=subInd

    % Outer grouping according to subjects
    n = n + 1;
    rowNames{n} = num2str(s);
    for m=1:length(fnames)    
        
        splits = cvll_basic(s).(fnames{m}).splits;
        nTrain = sum(splits(:,1),1);
        nTest = size(splits(:,1),1) - nTrain;
        
        % Only if the model is available
        if ~isempty(getfield(cvll_basic(s),fnames{m}))      
            
            % Get model evidence for given model
            cv_dist = getfield(getfield(cvll_basic(s),fnames{m}),'LLH');

            % Model evidence (test set)
            cv_perf(s).(fnames{m}) = median(cv_dist(:,2));             
            cv_SE(s).(fnames{m}) = iqr(bootstrp(100,@median,cv_dist(:,2)));     % error of the median

            % Overfitting
            % Median difference between per trial LL
            % Multiplied by number of test set
            cv_overfit(s).(fnames{m}) = median(cv_dist(:,1)./nTrain-cv_dist(:,2)./nTest)*nTest;

            % R^2
            tmp = getfield(getfield(cvll_basic(s),fnames{m}),'Rsq');
            cv_Rsq(s).(fnames{m}) = median(tmp(:,2));

            clear res cv_dist tmp train test
        else
            % Fill up with NaN
            cv_perf(s).(fnames{m}) = nan;
            cv_SE(s).(fnames{m}) = nan;
            cv_overfit(s).(fnames{m}) = nan;
            cv_Rsq(s).(fnames{m}) = nan;
        end
    end
end
clear n s m

%% Make arrays and tables
% Model evidence
M = squeeze(cell2mat(struct2cell(cv_perf)))';
tM = array2table(M,'RowNames',rowNames,'VariableNames',fnames);

% Standard Error of model evidence
SE = squeeze(cell2mat(struct2cell(cv_SE)))';
tSE = array2table(SE,'RowNames',rowNames,'VariableNames',fnames);

% R^2
Rsq = squeeze(cell2mat(struct2cell(cv_Rsq)))';
tRsq = array2table(Rsq,'RowNames',rowNames,'VariableNames',fnames);

% Index for over fitting
OF = squeeze(cell2mat(struct2cell(cv_overfit)))';
tOF = array2table(OF,'RowNames',rowNames,'VariableNames',fnames);

%% Save
save('model_evidence', 'M', 'tM', 'tSE', 'tRsq', 'tOF', 'fnames', 'rowNames');