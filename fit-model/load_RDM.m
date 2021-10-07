function [out_data, out_tbl] = load_RDM(exp)


if nargin == 0
    exp = 'exp2';
end
exp


% data to load
switch exp
    
    
    case {'exp1', 'exp1_dyn'}
        
        
        % === load experiment 1
        
        pts     = [3:52, 1:9];
        vers    = [repmat({'dots2'}, 1, length(3:52)), repmat({'dots2a'}, 1,length(1:9))];
        
        
    case {'exp2', 'exp2_dyn'}
        
        
        % === load experiment 2
        
        pts     = setxor(2:89, 32);
        vers    = repmat({'colmo'}, 1, length(pts));
        
    case {'exp3', 'exp3_dyn'}
        
        
        % === load experiment 2
        
        pts = setxor([100:200], 105);
        vers    = repmat({'colmo'}, 1, length(pts));
        
    case {'exp23', 'exp23_dyn'}
        
        
        % === load experiment 2-3
        
        pts     = [setxor(2:89, 32), setxor([100:200], 105)];
        vers    = repmat({'colmo'}, 1, length(pts));
        
end


resultsDir = '../data' % results directory


%% analysis parameters

plotClean   = 0
shiftRT     = 0.20
ptMinAcc    = 0.70


%% extract data & make table
[rt0, rt1, crt0, crt1, cslrt0, cslrt1, ...
    vrt0, vrt1, srt0, srt1,...
    acc0, acc1, ...
    resp0, resp1, respR, ...
    color, motion,...
    dist0, dist1, ...
    d5Dist0, d5Dist1, d11Dist0, d11Dist1,nDist,...
    targ0, targ1, ...
    d5Targ0, d5Targ1, d11Targ0, d11Targ1,...
    pes, ...
    blockTrials, totalTrials, discTrials, nTrials, ...
    splitHalf, rewBlock, blockNum, ...
    colBlkInd, ordType, expType,...
    pt] = deal([]);

maxRT = 0; ITI = [0,0];

[pt_age, pt_sex, pt_bonus] = deal([]);

vec = @(x) x(:);
center = @(x) x - nanmean(x);

pc = 1;
for pp  = 1:length(pts)
    
    % load file
    switch vers{pp}
        case 'dots2'
            fn = dir([resultsDir, '/pt' num2str(pts(pp)), '*dots_2.mat']);
        case 'dots2a'
            fn = dir([resultsDir, '/*dots_2a*' num2str(pts(pp)), '.mat']);
        case 'colmo'
            fn = dir([resultsDir, '/*dots' num2str(pts(pp)), '.mat']) ;
    end
    
    if isempty(fn)
        continue;
    end
    
    r = load(fullfile(resultsDir, fn(end).name));
    
    
    % demographics
    pt_age = [pt_age; r.data.ptInfo{2}];
    pt_sex = [pt_sex; strcmp(r.data.ptInfo{3}, 'F')];
    if  strcmp(exp, 'exp3')
        pt_bonus = [pt_bonus; ceil(r.data.rew.totDol)];
    end
    
    % === loop over blocks ===================================================
    blkNames = intersect(fieldnames(r.data.results), {'train', 'main', 'rand', 'ord'});
    
    for bb = 1:length(blkNames)
        
        
        % organize data
        switch vers{pp}
            case 'colmo'
                sesh = r.data.results.(blkNames{bb}).session;
            otherwise
                sesh = r.data.results.(blkNames{bb}).behav;
        end
        
        
        % attend-color / attend-motion
        nTrials = [nTrials; unique(sum(isfinite(sesh.task.color)))];
        colBlk = sum(isfinite(sesh.task.color)) > 66;
        [blkSel0,blkSel1] = deal(isfinite(sesh.task.color));
        
        blkSel0(1,:) = 0;
        
        blkSel1(end,:) = 0;
        blkSel1(nanmax(sum(isfinite(sesh.task.color(:,~colBlk)))), ~colBlk) = 0;
        
        
        colBlkMx = repmat(sum(isfinite(sesh.behav.rt)) > 66, [size(sesh.behav.rt,1),1]);
        colBlkInd = [colBlkInd; colBlkMx(blkSel0)];
        
        
        % === RT
        %         p_rt0 = sesh.behav.rt(blkSel0);
        %         p_rt1 = sesh.behav.rt(blkSel1);
        %         trialZeros = zeros(length(p_rt0),1);
        
        trialZeros = zeros(length(sesh.behav.rt(blkSel0)),1);
        
        %         rt0 = [rt0; p_rt0];
        %         rt1 = [rt1; p_rt1];
        
        rt0 = [rt0; sesh.behav.rt(blkSel0)];
        rt1 = [rt1; sesh.behav.rt(blkSel1)];
        
        
        % should center seperate for each analysis :(
        %         cp_rt0 = p_rt0;
        %         cp_rt1 = p_rt1;
        %
        %         cp_rt0(cp_rt0 <= minRT | cp_rt0 >= maxRT | cp_rt1 <= minRT | cp_rt1 >= maxRT | colBlkInd == 0)  = nan; %only good trials w/o
        %         cp_rt1(cp_rt1 <= minRT | cp_rt1 >= maxRT | cp_rt1 <= minRT | cp_rt1 >= maxRT | colBlkInd == 0)  = nan;
        %
        %         % center RT
        %         crt0 = [crt0; center(cp_rt0)];
        %         crt1 = [crt1; center(cp_rt1)];
        %
        %         % center shifted log RT
        %         cslrt0 = [cslrt0; center(log(cp_rt0-shiftRT))];
        %         cslrt1 = [cslrt1; center(log(cp_rt1-shiftRT))];
        
        
        % === accuracy
        sesh.behav.acc(~isfinite(sesh.behav.acc)) = 0;
        acc0 = [acc0; (sesh.behav.acc(blkSel0))];
        acc1 = [acc1; (sesh.behav.acc(blkSel1))];
        
        
        % === reponse hand
        resp0 = [resp0; (sesh.task.corrResp(blkSel0))];
        resp1 = [resp1; (sesh.task.corrResp(blkSel1))];
        
        % reponse repeat/switch
        try
            respR = [respR; (sesh.task.corrResp(blkSel0)) == (sesh.task.corrResp(blkSel1))];
        catch
            keyboard
        end
        
        % === color
        switch vers{pp}
            case 'colmo'
                colorNames = r.data.results.rand.param.colorNames;
            otherwise
                colorNames = r.data.results.train.param.colorNames;
        end
        
        colNum = (sesh.task.color(blkSel0));
        color = [color; colorNames(colNum)'];
        
        
        % === motion
        motion = [motion; (sesh.task.motion(blkSel0) == 0)];
        
        
        % === distractor
        dlvl = [];
        if isfield(r.data.param, 'driftConfLevels') && strcmp(blkNames{bb}, 'main')
            distlist = linspace(-max(r.data.param.cohRange), max(r.data.param.cohRange), r.data.param.driftConfLevels)./1000;
            normDist = distlist(sesh.confs);
            dlvl = r.data.param.driftConfLevels;
        else
            distlist = linspace(-max(r.data.param.cohRange), max(r.data.param.cohRange), r.data.param.confLevels)./1000;
            normDist = distlist(sesh.confs);
            dlvl = [dlvl, r.data.param.confLevels];
        end
        
        if bb ==1
            nDist = [nDist; min(dlvl)];
        end
        
        
        dist0 = [dist0; -(normDist(blkSel0))];
        dist1 = [dist1; -(normDist(blkSel1))];
        
        d5Dist0 = [d5Dist0; discretize((normDist(blkSel0)), 1.1*linspace(min(distlist), max(distlist), 6))];
        d5Dist1 = [d5Dist1; discretize((normDist(blkSel1)), 1.1*linspace(min(distlist), max(distlist), 6))];
        
        d11Dist0 = [d11Dist0; discretize((normDist(blkSel0)), 1.1*linspace(min(distlist), max(distlist), 12))];
        d11Dist1 = [d11Dist1; discretize((normDist(blkSel1)), 1.1*linspace(min(distlist), max(distlist), 12))];
        
        
        %         distlist = linspace(-max(r.data.param.cohRange), max(r.data.param.cohRange), 15)./1000;
        %         da = discretize((distlist), 1.1*linspace(min(distlist), max(distlist), 12))
        
        
        % === targets
        switch vers{pp}
            
            case 'colmo'
                
                targlist = center(linspace(min(r.data.param.attendCoh), max(r.data.param.attendCoh), r.data.param.confLevels)./1000);
                normTarg = targlist(sesh.attendCohs);
                
                
                targ0 = [targ0; (normTarg(blkSel0))];
                targ1 = [targ1; (normTarg(blkSel1))];
                
                %                 d5Targ0 = [d5Targ0; discretize(tiedrank(normTarg(blkSel0)), 5)];
                %                 d5Targ1 = [d5Targ1; discretize(tiedrank(normTarg(blkSel1)), 5)];
                %
                %                 d11Targ0 = [d11Targ0; discretize(tiedrank(normTarg(blkSel0)), 11)];
                %                 d11Targ1 = [d11Targ1; discretize(tiedrank(normTarg(blkSel1)), 11)];
                
                d5Targ0 = [d5Targ0; discretize((normTarg(blkSel0)), 1.1*linspace(min((targlist)), max((targlist)), 6))];
                d5Targ1 = [d5Targ1; discretize((normTarg(blkSel1)), 1.1*linspace(min((targlist)), max((targlist)), 6))];
                
                d11Targ0 = [d11Targ0; discretize((normTarg(blkSel0)), 1.1*linspace(min((targlist)), max((targlist)), 12))];
                d11Targ1 = [d11Targ1; discretize((normTarg(blkSel1)), 1.1*linspace(min((targlist)), max((targlist)), 12))];
                
                
                %                 targlist = (linspace(min(r.data.param.attendCoh), max(r.data.param.attendCoh), r.data.param.confLevels)./1000);
                %                 ta = discretize(center(targlist), 1.1*linspace(min(center(targlist)), max(center(targlist)), 6))

                
                
            otherwise
                
                targ0 = [targ0; nan*(normDist(blkSel0))];
                targ1 = [targ1; nan*(normDist(blkSel1))];
                
                d5Targ0 = [d5Targ0; nan*(normDist(blkSel0))];
                d5Targ1 = [d5Targ1; nan*(normDist(blkSel1))];
                
                d11Targ0 = [d11Targ0; nan*(normDist(blkSel0))];
                d11Targ1 = [d11Targ1; nan*(normDist(blkSel1))];
                
        end
        
        
        % === PES
        pesAcc = [sesh.behav.acc; nan(1, size(blkSel0,2))];
        
        
        blkSelPES0 = blkSel0;
        blkSelPES0(2,:) = 0;
        
        blkSelPES0 = [blkSelPES0; true(1, size(blkSel0,2))];

        
        blkSelPES1 = blkSel1;
        blkSelPES1(max(sum(isfinite(sesh.task.color(:,~colBlk))))-1, ~colBlk) = 0;
        blkSelPES1(max(sum(isfinite(sesh.task.color(:,colBlk))))-1, colBlk) = 0;
        
        blkSelPES1 = [blkSelPES1; true(1, size(blkSel0,2))];

        
        
        pes = [pes; pesAcc(blkSelPES0) - pesAcc(blkSelPES1)];
        
        
        % === trial count
        ptBlockTrials = cumsum(isfinite(sesh.task.color));
        normBlockTrials = (ptBlockTrials ./ max(ptBlockTrials)) - .5;
        ptTotalTrials = [1:length(trialZeros)]';
        normTotalTrials = (ptTotalTrials ./ max(ptTotalTrials)) - .5;
        
        blockTrials = [blockTrials; normBlockTrials(blkSel0)];
        totalTrials = [totalTrials; normTotalTrials];
        discTrials = [discTrials; discretize(tiedrank(normTotalTrials), 11) - 6];
        
        
        % === block type
        blockCount = ones(size(blkSel0));
        blockCount = blockCount .* [1:size(blockCount,2)];
        blockNum   = [blockNum;  center(blockCount(blkSel0))];
        
        
        % === reward
        try
            rew = zeros(size(blkSel0));
            rew(:, setxor(1:24,r.data.param.rewBlock)) = 1;
            rew(:, r.data.param.rewBlock) = 2;
            rewBlock = [rewBlock; rew(blkSel0)];
            
        catch
            
            rewBlock = [rewBlock; trialZeros*nan];
            
        end
        
        
        % === experiment
        switch exp
            case {'exp23', 'exp23_dyn'}
                if pts(pp) < 100
                    expType = [expType; trialZeros + 2];
                else
                    expType = [expType; trialZeros + 3];
                end
                
            otherwise
                expType = [expType; trialZeros + str2double(exp(4))];
        end
        
        
        % === session type
        switch vers{pp}
            case 'colmo'
                ordType = [ordType; repmat(blkNames(bb), length(trialZeros),1)];
            otherwise
                switch blkNames{bb}
                    case 'train'
                        ordType = [ordType; repmat({'rand'}, length(trialZeros),1)];
                    case 'main'
                        ordType = [ordType; repmat({'ord'}, length(trialZeros),1)];
                end
        end
        
        % === participant
        pt = [pt; trialZeros + pc];
        
    end
    
    pc = pc+1;
    
    
    maxRT(pp) = r.data.param.CTS.MaxRT;
    ITI(pp,:) = r.data.param.ITI;

end

exp
ITI
maxRT
mean(nTrials) 
range(nTrials)
% keyboard;


%% make table


tbl = table;

tbl.pt          = pt;
tbl.expType     = expType;
tbl.ordType     = ordType;
tbl.colBlk      = colBlkInd;

tbl.rt0         = rt0;
tbl.rt1         = rt1;
% tbl.crt0        = crt0;
% tbl.crt1        = crt1;
% tbl.cslrt0      = cslrt0;
% tbl.cslrt1      = cslrt1;
% tbl.vrt0        = vrt0;
% tbl.vrt1        = vrt1;
% tbl.srt0        = srt0;
% tbl.srt1        = srt1;

tbl.acc0        = acc0;
tbl.acc1        = acc1;

tbl.resp0       = resp0;
tbl.resp1       = resp1;
tbl.respR       = respR;

tbl.color       = color;
tbl.motion      = motion;

tbl.dist0       = dist0;
tbl.dist1       = dist1;
tbl.d5Dist0     = d5Dist0;
tbl.d5Dist1     = d5Dist1;
tbl.d11Dist0    = d11Dist0;
tbl.d11Dist1    = d11Dist1;

tbl.targ0       = targ0;
tbl.targ1       = targ1;
tbl.d5Targ0     = d5Targ0;
tbl.d5Targ1     = d5Targ1;
tbl.d11Targ0    = d11Targ0;
tbl.d11Targ1    = d11Targ1;

tbl.pes         = pes;

tbl.blockTrials = blockTrials;
tbl.totalTrials = totalTrials;
tbl.discTrials  = discTrials;
tbl.blockNum    = blockNum;

tbl.rew         = rewBlock;



%% make additional variables

% log rt
tbl.lrt0        = log(rt0);
tbl.lrt1        = log(rt1);


% shifted log rt
tbl.slrt0        = log(rt0-shiftRT);
tbl.slrt1        = log(rt1-shiftRT);


% error
tbl.err0        = ~acc0;
tbl.err1        = ~acc1;


% modified distractor
tbl.binDist0    = sign(dist0);
tbl.absDist0    = abs(dist0);

tbl.posDist0    = dist0 .* (dist0 > 0);
tbl.negDist0    = dist0 .* (dist0 < 0);

tbl.binDist1    = sign(dist1);
tbl.absDist1    = abs(dist1);

tbl.posDist1    = dist1 .* (dist1 > 0);
tbl.negDist1    = dist1 .* (dist1 < 0);


% target
tbl.binTarg0    = sign(targ0);
tbl.binTarg1    = sign(targ1);


% log target
tbl.ltarg0      = log(abs(.6 + targ0)) - mean(unique(log(abs(.6 + targ0))));
tbl.ltarg1      = log(abs(.6 + targ1)) - mean(unique(log(abs(.6 + targ1))));



% response coding
tbl.corrResp = tbl.resp0;
tbl.corrResp(tbl.corrResp == 1) = 0;
tbl.corrResp(tbl.corrResp == 3) = 1;

tbl.sgnResp = tbl.corrResp;
tbl.sgnResp(tbl.sgnResp==0) = -1;

tbl.ptResp = tbl.corrResp;
tbl.ptResp(tbl.acc0 == 0) = tbl.ptResp(tbl.acc0 == 0) == 0;


tbl.corrResp1 = tbl.resp1;
tbl.corrResp1(tbl.corrResp1 == 1) = 0;
tbl.corrResp1(tbl.corrResp1 == 3) = 1;

tbl.sgnResp1 = tbl.corrResp1;
tbl.sgnResp1(tbl.sgnResp1==0) = -1;

tbl.ptResp1 = tbl.corrResp1;
tbl.ptResp1(tbl.acc1 == 0) = tbl.ptResp1(tbl.acc1 == 0) == 0;

tbl.sgnPtResp1 = tbl.ptResp1;
tbl.sgnPtResp1(tbl.sgnPtResp1==0) = -1;

%% categorical variables
tbl.resp0 = categorical(tbl.resp0, [1,3], {'L', 'R'});
tbl.resp1 = categorical(tbl.resp1, [1,3], {'L', 'R'});
tbl.respR = categorical(tbl.respR, [0,1], {'switch', 'repeat'});

if any(isnan(tbl.rew))
    tbl.rew = categorical(tbl.rew, [nan,1,2], {'low', 'low', 'high'});
else
    tbl.rew = categorical(tbl.rew, [1,2], {'low', 'high'});
end

% tbl.ordType = categorical(tbl.ordType, [1, 2], {'rand', 'pred'});
tbl.ordType = categorical(tbl.ordType);

tbl.colBlk = categorical(tbl.colBlk, [0, 1], {'motion', 'color'});


tbl.color = categorical(tbl.color);
tbl.motion = categorical(tbl.motion, [0,1], {'left', 'right'});


switch exp
    case {'exp23', 'exp23_dyn'}
        tbl.expType = categorical(tbl.expType, [2, 3], {'2', '3'});
end


%% clean data

%% if any(tbl.expType == 3)
%
%     grpTbl = grpstats(tbl, {'pt', 'rew'}, {'mean', 'median'}, 'DataVars', {'acc0', 'rt0'});
%     grpTbl = grpstats(grpTbl, {'pt'}, {'min'}, 'DataVars', {'mean_acc0', 'median_rt0'});
%
%     grpTbl.mean_acc0  = grpTbl.min_mean_acc0;
%     grpTbl.median_rt0 = grpTbl.min_median_rt0;
%
% elseif any(tbl.expType == 1)
%
%     grpTbl = grpstats(tbl, {'pt', 'ordType'}, {'mean', 'median'}, 'DataVars', {'acc0', 'rt0'});
%     grpTbl = grpstats(grpTbl, {'pt'}, {'min'}, 'DataVars', {'mean_acc0', 'median_rt0'});
%
%     grpTbl.mean_acc0  = grpTbl.min_mean_acc0;
%     grpTbl.median_rt0 = grpTbl.min_median_rt0;
%
% else
%
%     grpTbl = grpstats(tbl, {'pt'}, {'mean', 'median'}, 'DataVars', {'acc0', 'rt0'});
%
% end

sel = ismember(tbl.colBlk, 'color') & ismember(tbl.ordType, 'rand');

grpTbl = grpstats(tbl(sel, :), {'pt'}, {'mean', 'median'}, 'DataVars', {'acc0', 'rt0'});


% plot grp stats
reportTbl = grpstats(tbl(sel, :), 'pt', 'mean', 'DataVars', 'acc0');

disp(reportTbl)

if plotClean == 1
    
    nexttile;
    histogram(grpTbl.mean_acc0, 25);
    xline(.75, '--k', 'LineWidth', 1);
    title(exp);
    xlim([.5 1]);
    
    keyboard;
    
end



% == select out good participants
% 1) mean accuracy > ptMinAcc
% 2) distractor levels > 2 (coding error)
% 3) completed half the experimeriment

goodSel = ...
    (grpTbl.mean_acc0 >  ptMinAcc) & ...
    (nDist > 2) & ...
    (reportTbl.GroupCount > median(reportTbl.GroupCount)/2);

goodPt  = grpTbl.pt(goodSel);
fprintf('good pt: %d \n', length(goodPt))
fprintf('removed pt: %d \n', length(grpTbl.mean_acc0) - length(goodPt))

tbl = tbl(ismember(tbl.pt, goodPt), :);


demo_tbl = table;
demo_tbl.age = pt_age(goodPt);
demo_tbl.sex = pt_sex(goodPt);
grpstats(demo_tbl, {}, {'mean', 'std', 'min', 'max'}, 'DataVars', {'age'})
grpstats(demo_tbl, {}, {'sum'}, 'DataVars', {'sex'})

if  strcmp(exp, 'exp3')
    demo_tbl.bonus = pt_bonus(goodPt);
    grpstats(demo_tbl, {}, {'mean', 'std'}, 'DataVars', {'bonus'})
end

% plot table
% summary(tbl);
% head(tbl);





%% make data structures

for pp = 1:length(goodPt)
    
    out_data(pp) = table2struct(tbl(ismember(tbl.pt, goodPt(pp)),:),'ToScalar',true);

end


out_tbl = tbl;
