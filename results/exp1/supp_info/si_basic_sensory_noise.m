% Supp Fig 5a

clear variables
close all

%% Setup
% Add path to auxiliary functions
addpath('.\..\..\..\src\');

% Subselect subjects
subInd = [1:24];

% Load data
load('.\..\..\..\data\exp1_data.mat');

setN = [3 5 7 8 9 10 11 12 13]'; % possible sample size

allNdiff = [];
allN = [];
alldecision = [];

%%
for s=subInd
    clear trials
    trials = trialData{s};
    
    [w,~, output(s)] = glmfit(trials.meanEvidence-0.5, trials.decision, 'binomial','link','probit');
    sub(s).ncdf_std = 1/w(2); % gaussian noise is inverse of weight
    sub(s).bias = w(1);
    
    % percentage of misclassified
    misclassified = 1*((trials.meanEvidence>0.5) ~= trials.decision);
    misclassified(trials.meanEvidence==0.5) = nan; % equal number of dots
    p_misclassified(s) = nanmean(misclassified);
    
    %
    N = trials.sampleSize;
    Nb = trials.meanEvidence.*N;
    Ndiff = Nb-(N-Nb); % difference of dots
    allN = [allN; N];
    allNdiff = [allNdiff;Ndiff];
    alldecision = [alldecision;trials.decision];
    
    
    ferr = @(q) negllh_pnoise(q, Ndiff, N, trials.decision);
opts = optimoptions('fmincon','Algorithm','interior-point','Display','none');
[wboth(:,s),~,EXITFLAG_both(s),output_both(s),~,~,H(:,:,s)] = fmincon(ferr,[1 1],[],[],[],[],[0 0],[],[],opts);


end

%% Aggregate across participants
Q = [sub.ncdf_std]';
sdMedian = median(Q);
bs = bootstrp(1000,@median,Q);
quantile(bs,[0.05 0.95]);

fprintf('- [result] estimated Gaussian SD = %.4f, 95%%-CI = (%.4f, %.4f)\n', median(Q), quantile(bs,[0.05 0.95]));

fprintf('- [result] misclassification rate: average = %.4f, standard deviation %.4f\n', mean(p_misclassified), std(p_misclassified));

%% Misclassification rate curve
for n=1:length(setN)  % for each sample length
    % Mask
    H = 0:setN(n); % possible number of heads (blue?)
    H(2*H==setN(n)) = []; % exclude situation with same number of dots of both colors
    T = setN(n) - H; % corresponding number of tails
    
    
    
    bin{n}.x = H./setN(n); % corresponding proportion of blue samples
    % N-specific bins
    for b=1:length(H) % for each proportion of blue samples
        R = []; % concatenate from all subjects
        
        for s=subInd % for each subject
            trials = trialData{s};
            
            mEv = trials.meanEvidence;
            N = trials.sampleSize;
            misclassified = (mEv>0.5) ~= trials.decision;
            
            assert(mean(unique(N)==setN)==1); % check that values used for number of dots do correspond
            
            % select trials from subject with given number of samples and
            % given proportion of blue dots
            nMask = N == setN(n);
            mEvMask = mEv == H./(H+T);
            mask = nMask & mEvMask;
            
            Rind = misclassified(mask(:,b)); % to plot single subject
            
            R = [R;Rind];
            
        end
        
        
        % average across subjects
        bin{n}.mean(b,1) = mean(R);
        bin{n}.SE(b,1) = std(R)/sqrt(length(R));
        
    end
    
    
end

%% Plot
width = 8;
height = 6;
LW = 1.2;
FS = 11;
figure(1);
figname = 'Misclassification_rate';
set(gcf,'name',figname);
clf;
hold on

colGrad = @(q) hsv2rgb([interp1([0 1]',[.08 .33]',q(:)), ones(numel(q),1)*1, ones(numel(q),1)*0.85]);
colGradModel = @(q) hsv2rgb([interp1([0 1]',[.08 .33]',q(:)), ones(numel(q),1)*.35, ones(numel(q),1)*1]);

line([1 1]*50,[0 1],'LineStyle','--','LineWidth',0.8,'Color',[1 1 1]*0.8);

for n=1:length(setN)
    h(n) = errorbar(bin{n}.x*100,100*bin{n}.mean,100*bin{n}.SE,'Color',colGrad((n-1)/(numel(setN)-1)),'LineWidth',LW,'LineStyle','-','CapSize',0);
end

xlim([0 100]);
%ylim([0.5 1]);

xlabel('proportion of blue samples (\%)', 'FontSize', FS, 'FontName', 'Times', 'Interpreter', 'latex');
ylabel('Misclassification rate (\%)', 'FontSize', FS, 'FontName', 'Times', 'Interpreter', 'latex');

set(gcf,'Color',[1,1,1]);

% Position plot on the screen for drawing
set(gcf, 'Units', 'centimeters', 'Position', [2 4 width height]);

% Position plot on the paper for printing
set(gcf, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'manual',...
    'PaperSize', [width height], 'PaperPosition', [0 0 width height]);

% Axes
set(gca, 'Box', 'off', 'FontSize', FS, 'FontName', 'Times', 'TickDir', 'out', 'OuterPosition', [0 0 1 1],...  % try to place axes first
    'XGrid', 'off',  'YGrid', 'off', 'Layer', 'top');

%% Print
filename = fullfile( '.\..\..\..\plots\exp1\supp_info', [figname '.png']);
print(gcf, '-dpng', '-r400', filename);


%% Fit model with concatenated data to split into sample size independent and sample size dependent


ferr = @(q) negllh_pnoise(q, allNdiff, allN, alldecision);
opts = optimoptions('fmincon','Algorithm','interior-point','Display','none');
[w_all,FVAL,EXITFLAG,output_all,~,~,H] = fmincon(ferr,[1 1],[],[],[],[],[0 0],[],[],opts);



function negLLH = negllh_pnoise(sigma, Ndiff, N, decision)
sig_ind = sigma(1); % sample size independent noise (std)
sig_dep = sigma(2); % sample size dependent noise (std)
noise_all = sqrt(sig_ind^2 + sig_dep^2*N.^2); % overall noise for each trial (std)
y = normcdf( Ndiff./ noise_all); % probability of selecting right response
LLH = sum(log(y(decision==1))) + sum(log(1-y(decision==0))); % loglikelihood
negLLH = -LLH;

end

