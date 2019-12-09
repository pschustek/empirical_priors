% Supp Fig 5b
clear variables
close all

%% Setup
% Add path to auxiliary functions
addpath('.\..\..\..\src\');

% Subselect subjects
subInd = [1:24];

% Load data
load('.\..\..\..\data\exp2_data.mat');

setN = [3:11]'; % possible sample size

allNdiff = [];
allN = [];
alldecision = [];

%%
for s=subInd
    clear trials
    trials = trialData{s};
    
    % percentage of misclassified
    evidenceopposing = 1*((trials.meanEvidence>0.5) ~= trials.decision);
    evidenceopposing(trials.meanEvidence==0.5) = nan; % equal number of dots
    p_misclassified(s) = nanmean(evidenceopposing);
    
    %
    N = trials.sampleSize;
    Nb = trials.meanEvidence.*N;
    Ndiff = Nb-(N-Nb); % difference of dots
    allN = [allN; N];
    allNdiff = [allNdiff;Ndiff];
    alldecision = [alldecision;trials.decision];
    
end


%% Evidence opposing rate curve
%for c=1:2
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
            evidenceopposing = (mEv>0.5) ~= trials.decision;
            
            assert(mean(unique(N)==setN)==1); % check that values used for number of dots do correspond
            
            % aligned mean evidence
            block = 2*(trials.blockBias-0.5);
            mEv = block.*(mEv-0.5) + 0.5;
            
            % select trials from subject with given number of samples and
            % given proportion of blue dots
            nMask = N == setN(n);
            mEvMask = mEv == H./(H+T);
            mask = nMask & mEvMask;
            
            
            Rind = evidenceopposing(mask(:,b)); % to plot single subject
            
            R = [R;Rind];
            
        end
        
        
        % average across subjects
        bin{n}.mean(b,1) = mean(R);
        bin{n}.SE(b,1) = std(R)/sqrt(length(R));        
    end
    
    
end
%end


%% Plot
width = 1*8;
height = 6;
LW = 1.2;
FS = 11;
figure(1);
figname = 'Evidence_opposing_choices';
set(gcf,'name',figname);
clf;

%for c=1:2
%subplot(1,2,c);
hold on

colGrad = @(q) hsv2rgb([interp1([0 1]',[.08 .33]',q(:)), ones(numel(q),1)*1, ones(numel(q),1)*0.85]);
colGradModel = @(q) hsv2rgb([interp1([0 1]',[.08 .33]',q(:)), ones(numel(q),1)*.35, ones(numel(q),1)*1]);

line([1 1]*50,[0 1],'LineStyle','--','LineWidth',0.8,'Color',[1 1 1]*0.8);

for n=1:length(setN)
       h(n) = errorbar(bin{n}.x*100,100*bin{n}.mean,100*bin{n}.SE,'Color',colGrad((n-1)/(numel(setN)-1)),'LineWidth',LW,'LineStyle','-','CapSize',0);
   % h(n) = errorbar(bin{n}.x*100,100*bin{n}.mean(:,c),100*bin{n}.SE(:,c),'Color',colGrad((n-1)/(numel(setN)-1)),'LineWidth',LW,'LineStyle','-','CapSize',0);
end

xlim([0 100]);
ylim([0 50]);

plot([50 50], ylim,'color',.7*[1 1 1]);

xlabel('Aligned stimulus evidence (\%)', 'FontSize', FS, 'FontName', 'Times', 'Interpreter', 'latex');
%if c==1
    ylabel('Evidence opposing choices (\%)', 'FontSize', FS, 'FontName', 'Times', 'Interpreter', 'latex');
%end
%end

set(gcf,'Color',[1,1,1]);

% Position plot on the screen for drawing
set(gcf, 'Units', 'centimeters', 'Position', [2 4 width height]);

% Position plot on the paper for printing
set(gcf, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'manual',...
    'PaperSize', [width height], 'PaperPosition', [0 0 width height]);

% Axes
set(gca, 'Box', 'off', 'FontSize', FS, 'FontName', 'Times', 'TickDir', 'out', ...  % try to place axes first
    'XGrid', 'off',  'YGrid', 'off', 'Layer', 'top');

%% Print
filename = fullfile( '.\..\..\..\plots\exp2\supp_info', [figname '.png']);
print(gcf, '-dpng', '-r400', filename);