%% FFIc simulation



%% init
clear; clc; close all;

addpath(genpath('../'));



% functions
center      = @(x) x-nanmean(x);

getBetaMean = @(exp, name) exp.results.alpha(ismember(exp.results.paramName, name));
getBetaSE   = @(exp, name) exp.results.stats.groupmeanerr(ismember(exp.results.paramName, name));
getBetaCI   = @(exp, name) exp.results.stats.groupmeanerr(ismember(exp.results.paramName, name)) .* abs(tinv(.025, size(exp.data,2)-exp.results.fitparams.Np));
getBetaAll  = @(exp, name) exp.results.link{ismember(exp.results.paramName, name)}(exp.results.E(ismember(exp.results.paramName, name),:))';
getBetaV    = @(exp, name) exp.results.V(ismember(exp.results.paramName, name),:)';
getBIC      = @(exp) exp.results.bf.bic';
getAIC      = @(exp) exp.results.bf.aic';





%% ====== LOAD DYNAMICS =========================================================
clear col polyMdl expMdl mot wn_col wn_mot




% tanh AR softplus
col{1} = load('2021-04-26_16-36_exp1_withinDualDyn.mat');
col{2} = load('2021-01-17_20-55_exp2_withinDualDyn.mat');
col{3} = load('2021-01-17_22-52_exp3_withinDualDyn.mat');



mot{1} = load('2021-03-08_05-53_exp1_withinMotionDualDyn.mat');
mot{2} = load('2021-03-07_23-36_exp2_withinMotionDualDyn.mat');
mot{3} = load('2021-03-08_05-59_exp3_withinMotionDualDyn.mat');


% tanh AR softplus
wn_col{1} = load('2021-01-17_21-31_exp1_withinStatic.mat');
wn_col{2} = load('2021-01-17_18-43_exp2_withinStatic.mat');
wn_col{3} = load('2021-01-20_00-43_exp3_withinStatic.mat');


% tanh AR softplus
wn_mot{1} = load('2021-01-17_21-40_exp1_motionStatic.mat');
wn_mot{2} = load('2021-01-17_22-28_exp2_motionStatic.mat');
wn_mot{3} = load('2021-01-19_10-12_exp3_motionStatic.mat');







disp('');
disp('loaded');
disp('');



% plot stats
disp_mdls = {'col', 'mot'}

% for cc = 1:3
%
%     npt(cc) = length(col{cc}.data);
%
%     for mm = 1:length(disp_mdls)
%
%         % == color
%         tbl = table;
%
%         tbl.pname   = eval([disp_mdls{mm},'{cc}.results.paramName'])';
%         tbl.df      = eval(['repmat(length(',disp_mdls{mm},'{cc}.data) - length(', disp_mdls{mm},'{cc}.results.stats.tval), [ length(', disp_mdls{mm},'{cc}.results.stats.tval), 1])']);
%         tbl.t       = eval([disp_mdls{mm},'{cc}.results.stats.tval']);
%         tbl.p       = eval([disp_mdls{mm},'{cc}.results.stats.p']);
%         tbl.d       = eval([disp_mdls{mm},'{cc}.results.stats.alpha./sqrt(diag(', disp_mdls{mm}, '{cc}.results.stats.groupvar))']);
%
%         fprintf('\n =============== exp %d %s =============== \n\n', cc, disp_mdls{mm});
%         disp(tbl)
%
%     end
%
%      fprintf('\n\n\n\n\n');
%
% end




% ==== combined stats
aaallT = struct;
for mm = 1:length(disp_mdls)
    
    % combined table
    allpname = [];
    for cc = 1:3
        allpname =  [allpname; eval([disp_mdls{mm},'{cc}.results.paramName'])'];
    end
    unqpname = unique(upper(allpname));
    alluname{mm} = unqpname;
    
    [dstat] = nan(length(unqpname), 3, 3);
    
    clear df
    
    
    for cc = 1:3
        

        pname   = eval([disp_mdls{mm},'{cc}.results.paramName'])';
        
        for uu = 1:length(unqpname)
            
            idx = ismember(upper(pname), unqpname{uu});
            
            if ~any(idx) | (any(regexpi(unqpname{uu}, 'targ')) & cc==1); continue; end
            
            
            dstat(uu,1,cc)      = eval([disp_mdls{mm},'{cc}.results.stats.tval(idx)']); % tval
            dstat(uu,2,cc)      = eval([disp_mdls{mm},'{cc}.results.stats.p(idx)']);    % pval
            dstat(uu,3,cc)      = eval([disp_mdls{mm},'{cc}.results.stats.alpha(idx)./sqrt(', disp_mdls{mm}, '{cc}.results.stats.groupvar(idx,idx))']); %cohen d
            
        end
        
        df(cc)      = eval(['length(',disp_mdls{mm},'{cc}.data) - length(', disp_mdls{mm},'{cc}.results.stats.tval)']);
        sqn(cc) = sqrt(eval(['length(',disp_mdls{mm},'{cc}.data)']));
        
    end
    
    
    dvals = squeeze(dstat(:,3,:));
    pvals = squeeze(dstat(:,2,:));
    sgnD = sign(dvals).*sign(dvals(:,end));
    
%     combP1 = 2*normcdf(nansum(norminv(1-((realmin+pvals)./2)).*df, 2) ./ sqrt(isfinite(pvals)*(df.^2)'), 'upper');
    combP = normcdf(nansum(sgnD.*-norminv(pvals).*sqn, 2) ./ sqrt(isfinite(pvals)*(sqn.^2)'), 'upper');
    
    diffSign = nanvar(sign(dvals),[],2) ~= 0;
    
    aaallT(mm).t = array2table([round(1e5*reshape(squeeze(dstat(:,3,:)), [length(unqpname), 3]))/1e5, combP, diffSign],...
        'VariableNames',...
        {'exp1 d',...
        'exp2 d', ...
        'exp3 d',...
        'combined p', 'diffSign'}, ...
        'RowNames', unqpname);
    
    disp(' ');
    disp(disp_mdls{mm});
    disp(aaallT(mm).t)
    disp(unqpname)
    disp(df)
    disp(' ');
    
    
end








% ====== compare exp 2 & 3
clear npt;
for ii = 1:3
            npt(cc) = length(col{cc}.data);

end

tcrt = ismember(col{2}.results.paramName, 'acc_{targCrt}');
cv2 = col{2}.results.stats.groupmeancovariance(tcrt,tcrt);
cv3 = col{3}.results.stats.groupmeancovariance(tcrt,tcrt);

targ_tval = (col{2}.results.alpha(tcrt) - col{3}.results.alpha(tcrt))./...
    sqrt(cv2 + cv3);

targ_df = (cv2 + cv3).^2 ./...
    ((cv2.^2)./(npt(2)-1) + (cv3.^2)./(npt(3)-1));

targ_pval = 2*tcdf(abs(targ_tval), targ_df, 'upper');





tcrt = ismember(col{2}.results.paramName, 'acc_{distCrt}');
cv2 = col{2}.results.stats.groupmeancovariance(tcrt,tcrt);
cv3 = col{3}.results.stats.groupmeancovariance(tcrt,tcrt);

dist_tval = (col{2}.results.alpha(tcrt) - col{3}.results.alpha(tcrt))./...
    sqrt(cv2 + cv3);

dist_df = (cv2 + cv3).^2 ./...
    ((cv2.^2)./(npt(2)-1) + (cv3.^2)./(npt(3)-1));

dist_pval = 2*tcdf(abs(dist_tval), dist_df, 'upper');


fprintf('\nexp2&3 targ dyn: t(%.5g)=%.5g, p=%.5g',  targ_df, targ_tval, targ_pval)
fprintf('\nexp2&3 dist dyn: t(%.5g)=%.5g, p=%.5g\n\n',  dist_df, dist_tval, dist_pval)








% ====== compare color & motion
% 
% for cc = 1:3
%     
%     npt = length(col{cc}.data);
%     
%     % acc targ
%     tcrt = ismember(col{cc}.results.paramName, 'acc_{targCrt}');
%     cv2 = col{cc}.results.stats.groupmeancovariance(tcrt,tcrt);
%     cv3 = mot{cc}.results.stats.groupmeancovariance(tcrt,tcrt);
%     
%     accTarg_tval = (col{cc}.results.alpha(tcrt) - mot{cc}.results.alpha(tcrt))./...
%         sqrt(cv2 + cv3);
%     
%     accTarg_df = (cv2 + cv3).^2 ./...
%         ((cv2.^2)./(npt-1) + (cv3.^2)./(npt-1));
%     
%     accTarg_pval = 2*tcdf(abs(accTarg_tval), accTarg_df, 'upper');
%     
%     
%     % acc targ
%     tcrt = ismember(col{cc}.results.paramName, 'acc_{distCrt}');
%     cv2 = col{cc}.results.stats.groupmeancovariance(tcrt,tcrt);
%     cv3 = mot{cc}.results.stats.groupmeancovariance(tcrt,tcrt);
%     
%     accDist_tval = (col{cc}.results.alpha(tcrt) - mot{cc}.results.alpha(tcrt))./...
%         sqrt(cv2 + cv3);
%     
%     accDist_df = (cv2 + cv3).^2 ./...
%         ((cv2.^2)./(npt-1) + (cv3.^2)./(npt-1));
%     
%     accDist_pval = 2*tcdf(abs(accDist_tval), accDist_df, 'upper');
%     
%     
%     
%     fprintf('\n exp %d', cc');
%     fprintf('\n == col-mot targCrt acc: t(%.5g)=%.5g, p=%.5g',  accTarg_df, accTarg_tval, accTarg_pval)
%     fprintf('\n == col-mot distCrt acc: t(%.5g)=%.5g, p=%.5g\n\n',  accDist_df, accDist_tval, accDist_pval)
%     
%     
% end





% ====== compare color & motion
clear pvals df rtDist_cvRatio sqrn sgnD
for cc = 1:3
    
    npt = length(col{cc}.data);
    cc_sqn = sqrt(npt);
    
    col_df  = npt - length(col{cc}.results.stats.tval);
    col_idf  = 1./col_df;
    mot_df  = npt - length(mot{cc}.results.stats.tval);
    mot_idf  = 1./mot_df;
    
    
    % acc targ
    tcrt = ismember(col{cc}.results.paramName, 'acc_{targCrt}');
    mtcrt = ismember(mot{cc}.results.paramName, 'acc_{targCrt}');
    
    cv2 = col{cc}.results.stats.groupmeancovariance(tcrt,tcrt);
    cv3 = mot{cc}.results.stats.groupmeancovariance(mtcrt,mtcrt);
    
    accTarg_tval = (col{cc}.results.alpha(tcrt) - mot{cc}.results.alpha(mtcrt))./...
        sqrt(cv2 + cv3);
    
    accTarg_df = (col_idf*cv2 + mot_idf*cv3).^2 ./...
        (col_idf.*(col_idf*cv2).^2 + mot_idf.*(mot_idf*cv3).^2);
    
    accTarg_pval = 2*tcdf(abs(accTarg_tval), accTarg_df, 'upper');
    
    
    % acc targ
    tcrt = ismember(col{cc}.results.paramName, 'acc_{distCrt}');
    mtcrt = ismember(mot{cc}.results.paramName, 'acc_{distCrt}');
    cv2 = col{cc}.results.stats.groupmeancovariance(tcrt,tcrt);
    cv3 = mot{cc}.results.stats.groupmeancovariance(mtcrt,mtcrt);
    if cc > 1
        accDist_cvRatio(cc) = (cv3/cv2);
    end
    accDist_tval = (col{cc}.results.alpha(tcrt) - mot{cc}.results.alpha(mtcrt))./...
        sqrt(cv2 + cv3);
    
    accDist_df = (col_idf*cv2 + mot_idf*cv3).^2 ./...
        (col_idf.*(col_idf*cv2).^2 + mot_idf.*(mot_idf*cv3).^2);
    
    accDist_pval = 2*tcdf(abs(accDist_tval), accDist_df, 'upper');
    
    
    
    % rt targ
    tcrt = ismember(col{cc}.results.paramName, 'rt_{targAcc}');
    mtcrt = ismember(mot{cc}.results.paramName, 'rt_{targAcc}');
    cv2 = col{cc}.results.stats.groupmeancovariance(tcrt,tcrt);
    cv3 = mot{cc}.results.stats.groupmeancovariance(mtcrt,mtcrt);
    
    if cc > 1
        rtTarg_cvRatio(cc) = (cv3/cv2);
    end
    
    rtTarg_tval = (col{cc}.results.alpha(tcrt) - mot{cc}.results.alpha(mtcrt))./...
        sqrt(cv2 + cv3);
    
    
    rtTarg_df = (col_idf*cv2 + mot_idf*cv3).^2 ./...
        (col_idf.*(col_idf*cv2).^2 + mot_idf.*(mot_idf*cv3).^2)
    
    rtTarg_pval = 2*tcdf(abs(rtTarg_tval), rtTarg_df, 'upper');
    
    
    % rt targ
    tcrt = ismember(col{cc}.results.paramName, 'rt_{distAcc}');
    mtcrt = ismember(mot{cc}.results.paramName, 'rt_{distAcc}');
    cv2 = col{cc}.results.stats.groupmeancovariance(tcrt,tcrt);
    cv3 = mot{cc}.results.stats.groupmeancovariance(mtcrt,mtcrt);
    
    if cc > 1
        rtDist_cvRatio(cc) = (cv3/cv2);
    end
    
    rtDist_tval = (col{cc}.results.alpha(tcrt) - mot{cc}.results.alpha(mtcrt))./...
        sqrt(cv2 + cv3);
    
    
    rtDist_df = (col_idf*cv2 + mot_idf*cv3).^2 ./...
        (col_idf.*(col_idf*cv2).^2 + mot_idf.*(mot_idf*cv3).^2)
    
    
    rtDist_pval = 2*tcdf(abs(rtDist_tval), rtDist_df, 'upper');
    
    
    fprintf('\n exp %d', cc');
    fprintf('\n == col-mot targ acc: t(%.5g)=%.5g, d=%.5g, p=%.5g',  accTarg_df, accTarg_tval, accTarg_tval./sqrt(accTarg_df), accTarg_pval)
    fprintf('\n == col-mot dist acc: t(%.5g)=%.5g, d=%.5g, p=%.5g',  accDist_df, accDist_tval, accDist_tval./sqrt(accDist_df), accDist_pval)
    fprintf('\n == col-mot targ rt: t(%.5g)=%.5g, d=%.5g, p=%.5g',  rtTarg_df, rtTarg_tval, rtTarg_tval./sqrt(rtTarg_df), rtTarg_pval)
    fprintf('\n == col-mot dist rt: t(%.5g)=%.5g, d=%.5g, p=%.5g\n\n',  rtDist_df, rtDist_tval, rtDist_tval./sqrt(rtDist_df), rtDist_pval)
    
    if cc==1
        [accTarg_tval, accTarg_pval, accTarg_df, rtTarg_tval, rtTarg_pval, rtTarg_df] =deal(nan);
    end
    
    
    pvals(:,cc) = [accTarg_pval;  accDist_pval; rtTarg_pval; rtDist_pval];
    df(:,cc) = [accTarg_df; accDist_df; rtTarg_df; rtDist_df];
    
    sqrn(:,cc) = [cc_sqn; cc_sqn; cc_sqn; cc_sqn];
    sgnD(:,cc) = sign([accTarg_tval; accDist_tval; rtTarg_tval; rtDist_tval]);
    
    
end

sgnD = sgnD == sign(sgnD(:,end));


% combP = 2*normcdf(nansum(norminv(1-((realmin+pvals)./2)).*df, 2) ./ sqrt(nansum((df.^2),2)), 'upper')
combP = normcdf(nansum(sgnD.*-norminv(pvals).*sqn, 2) ./ sqrt(isfinite(pvals)*(sqn.^2)'), 'upper')





%% ~~~~ Updated FFIc STATIC ~~~~~~~

nReals = 10000;          % how many simulated trials?



% ===== simulation parameters  =====
threshType = 'race';    % which threshold {race, diff}
doRelu = 1;             % force positive accumulators?

nChoices = 2;           % how many accumulators?
nDt    = 2000;          % max number of timesteps
dt = .001;              % simulation dt

% ===== task parameters =====

nTarg = 11;
nDist = 11;
targRange = linspace(.6,.95,nTarg);
distRange = linspace(-1,1,nDist);

[targ, dist] = meshgrid(targRange, distRange);
targ = tanh(1.5*targ(:))./tanh(1.5);
ntarg = tanh(1.5*(1-targ(:)))./tanh(1.5);
dist = tanh(1.35*dist(:))./tanh(1.35);

nConds = length(targ);





% ========= run COLOR STATIC sim




% decision dynamics

decay = 0;          % accumulator-specific decay
cross = 0;          % accumulator competition
ffinh = 1;          % feedforward inhibition
sd    = 1;          % diffusion noise
lapse = 0;          % lapse rate


A = cross*ones(nChoices) + (decay - cross)*eye(nChoices)
B = [...
    1, -ffinh;...
    -ffinh, 1;...
    ]



% gain dynamics
% targ0       = 1;
% targAim     = 6;
% targGain    = 3;
%  
% 
% dist0       = 1.25;
% distAim     = 0;
% distGain    = 4.5;
% 
% salSD       = 0.1;
% 
%  % non-decision time
% t0 = .300;         
% 
% % decision threshold
% thresh = max(.01, linspace(1.5, .1, nDt))';         






% gain dynamics
targ0       = 3.5;
dist0       = 1;

salSD       = 0;

 % non-decision time
t0 = .300;         

% decision threshold
thresh = max(.01, linspace(1.5, .1, nDt))';         
% thresh = max(.01, linspace(1.25, 1.25, nDt))';         
% thresh = max(.01, linspace(1, 1, nDt))';         




% a =  ones(2)*-25
% B = [...
%     1, 0;...
%     0, 1;...
%     ]




sddt = sd*sqrt(dt);
salSDdt = salSD*sqrt(dt);


% run sims  ======================

clear meanRT Accuracy static_allRT static_allChoice

tic
parfor cc = 1:nConds
    
    x = zeros(nChoices,nReals);
    all_x = nan(nDt,nChoices,nReals);
    
    %     targ0       = 2;
    %     targAim     = 5;
    %     targGain    = 1.5;
    %
    %     dist0       = 1;
    %     distAim     = 0;
    %     distGain    = 1.5;
    
    targSal = targ0;
    distSal = dist0;
    

    I = [...
        targ(cc),  abs(dist(cc))*(dist(cc)>0);...
        ntarg(cc), abs(dist(cc))*(dist(cc)<0);...
        ];
    
    
    for tt = 1:nDt
        
        % update gain
        %         targSal = targSal + targGain*(targAim - targSal)*dt + randn*salSDdt;
        %         distSal = distSal + distGain*(distAim - distSal)*dt + randn*salSDdt;
        
        % update model
        all_x(tt,:,:) = x;
        x = x + A*(x*dt) + B*(I*[targSal;distSal]*dt + randn(size(x))*sddt);
        
        % threshold
        if doRelu
            x = max(x,0);
        end
        
        % break if all accumulators reach threshold
        switch threshType
            
            case 'race'
                
                if ~any(max(x)<thresh(tt))
                    break;
                end
                
            case 'diff'
                
                if ~any(abs(diff(x))<thresh(tt))
                    break;
                end
                
        end
        
    end
    
    [RT, Choice] = deal(nan(nReals, 1));
    
    for rr = 1:nReals
        try
            switch threshType
                case 'race'
                    RT(rr)      = find(max(all_x(:,:,rr),[],2) >= thresh, 1);
                case 'diff'
                    RT(rr)      = find(abs(diff(all_x(:,:,rr),[],2)) >= thresh, 1);
            end
            
            [~,Choice(rr)]  = max(all_x(RT(rr),:,rr));
            
        catch
            RT(rr) = nan;
            Choice(rr) = 0;
        end
    end
    
    stdRT(cc) = nanstd(RT(Choice == 1)*dt + t0);
    meanRT(cc) = nanmean(RT(Choice == 1)*dt + t0);
    Accuracy(cc) = (1-lapse)*nanmean(Choice == 1) + lapse*.5;
    
    static_allRT{cc} = RT;
    static_allChoice{cc} = Choice;
    
end
toc


% get mean
rStdRTStatic = reshape(stdRT,[nTarg, nDist]);
rRTStatic = reshape(meanRT,[nTarg, nDist]);
rAccStatic = reshape(Accuracy,[nTarg, nDist]);




% make table
[t_dist, t_targ, t_acc, t_rt] = deal([]);
for cc = 1:nConds
    
    t_acc  = [t_acc; static_allChoice{cc}==1];
    t_rt   = [t_rt; static_allRT{cc}*dt + t0];
    t_targ = [t_targ; ones(nReals,1)*targ(cc)];
    t_dist = [t_dist; ones(nReals,1)*dist(cc)];
    
end



static_tbl = table;
static_tbl.acc = t_acc;
static_tbl.rt = t_rt;
static_tbl.lrt = log(t_rt);
static_tbl.crt = t_rt - nanmean(t_rt);
static_tbl.rt5 = discretize(tiedrank(t_rt),5);
static_tbl.targ = t_targ;
static_tbl.targ7 = discretize(tiedrank(t_targ),7);
static_tbl.dist = t_dist;
static_tbl.dist7 = discretize(tiedrank(t_dist),7);



% regression


% rt_mdl = fitlm(static_tbl, 'lrt ~ targ*dist', 'Exclude', static_tbl.acc == 0)
% acc_mdl = fitglm(static_tbl, 'acc ~ targ*dist', 'Distribution', 'binomial')
dyn_mdl = fitglm(static_tbl, 'acc ~ (targ + dist)*crt', 'Distribution', 'binomial')
dyn_mdl = fitlm(static_tbl, 'crt ~ (targ + dist)*acc')











 % PLOT DYN  



 







mdl = col;
% mdl = exp23;


[rt, rawRt, crt, d5rt, rt_pred, rt_targ, rt_dist, rt_acc, rt_pt,...
    acc, cacc, acc_pred, acc_targ, acc_dist, acc_rt, acc_pt,...
    all_ndt, rt_targ_all, rt_dist_all, acc_targ_all, acc_dist_all, ...
    ] = deal([]);

prevCor = [];

pt_count = 1;

nCohDisc = 7;


for mm = 1:3
    
    
    data = mdl{mm}.data;
    results = mdl{mm}.results;
    
    param = results.E';
    ndt = ...
        results.link{ismember(results.paramName, 'rt_{ndt}')} ...
        (param(:, ismember(results.paramName, 'rt_{ndt}')));
    
    all_ndt = [all_ndt; ndt];
    
    for tt = 1:size(data,2)
        
        %                 if mm == 1 && tt == size(data,2)
        %                     continue
        %                 end
        
        % simulation parameters
        [~, out]=eval([results.likfun, '(transpose(param(tt,:)), data(tt), [], [], 0, [])']);
        
        
        % rt
        prt         = log(out.rt.rt(out.rt.rt_sel) - ndt(tt));
        prt(out.rt.rt(out.rt.rt_sel) <= ndt(tt)) = nan;
        rt          = [rt; prt];
        crt         = [crt; prt - nanmean(prt)];
        rawRt       = [rawRt; out.rt.rt(out.rt.rt_sel)];

        
        rt_pred     = [rt_pred; out.rt.rt_yhat];
        
        rt_targ     = [acc_targ; discretize(tiedrank(out.data.targ0(out.rt.rt_sel)),nCohDisc)];
        rt_dist     = [acc_dist; discretize(tiedrank(out.data.dist0(out.rt.rt_sel)),nCohDisc)];
        
        
        rt_targ_all     = [rt_targ_all; out.data.d11Targ0(out.rt.rt_sel)];
        rt_dist_all     = [rt_dist_all; 12-out.data.d11Dist0(out.rt.rt_sel)];
        
        
        rt_acc      = [rt_acc; out.data.acc0(out.rt.rt_sel)];
        
        rt_pt       = [rt_pt; out.rt.rt_X(:,1).*pt_count];
        
        
        % acc
        d5rt          = [d5rt; discretize(tiedrank(out.data.rt0(out.acc.acc_sel)),5)];
        
        acc         = [acc; out.data.acc0(out.acc.acc_sel)];
        cacc        = [cacc; center(out.data.acc0(out.acc.acc_sel))];
        
        yhat        = [out.acc.acc_yhat];
        yhat(out.data.corrResp(out.acc.acc_sel)==0) = 1-yhat(out.data.corrResp(out.acc.acc_sel)==0);
        acc_pred    = [acc_pred; yhat];
        
        
        acc_targ     = [acc_targ; discretize(tiedrank(out.data.targ0(out.acc.acc_sel)),nCohDisc)];
        acc_dist     = [acc_dist; discretize(tiedrank(out.data.dist0(out.acc.acc_sel)),nCohDisc)];
        
        acc_targ_all     = [acc_targ_all; out.data.d11Targ0(out.acc.acc_sel)];
        acc_dist_all     = [acc_dist_all; 12-out.data.d11Dist0(out.acc.acc_sel)];
        
        acc_rt          = [acc_rt; out.data.rt0(out.acc.acc_sel)];
        acc_pt      = [acc_pt; out.acc.acc_X(:,1).*pt_count];
        
        pt_count = pt_count+1;
        
    end
    
    
end



rt_tbl = array2table([rt, rawRt, crt, rt_pred, rt_targ, rt_dist, rt_targ_all, rt_dist_all, rt_acc, rt_pt ], 'VariableNames', ...
    {'rt', 'rawRt', 'crt', 'rtY', 'targ', 'dist', 'targAll', 'distAll', 'acc', 'pt'});

acc_tbl = array2table([acc, cacc, acc_pred, acc_rt, d5rt, acc_targ, acc_dist,  acc_targ_all, acc_dist_all, acc_pt], 'VariableNames', ...
    {'acc', 'cacc', 'accY', 'rt', 'd5rt', 'targ', 'dist', 'targAll', 'distAll', 'pt'});







% ===== start plot
figure('units','inch','position',[0,0,1.5*7.5,6]); hold on;
tiledlayout(2,4, 'TileSpacing', 'compact', 'Padding', 'compact');


 
 

% ======= DYNAMICS





load('batlow.mat');
cols = batlow;
cIdx = round(linspace(1, length(cols), 7));

rtacc_cols = {[36, 123, 160]/255, [242, 95, 92]/255};





% ============== plot targ dynamics acc ==============
grpTbl = grpstats(acc_tbl, {'pt', 'd5rt', 'targ'});
grpTbl = grpstats(grpTbl, {'d5rt', 'targ'}, {'mean','sem'});


figP = nexttile([1,1]); hold on;

xoff = linspace(-.2,.2,5);
allY = [];
for ii = 1:5
    
    mAcc = mean(grpTbl.mean_mean_acc(grpTbl.d5rt == ii));
    
    e=errorbar([1:nCohDisc]+xoff(ii), (grpTbl.mean_mean_acc(grpTbl.d5rt == ii)), grpTbl.sem_mean_cacc(grpTbl.d5rt == ii),...
        'o', 'LineWidth', 1);
    e.Color = cols(cIdx(ii+1),:);
    e.MarkerFaceColor = cols(cIdx(ii+1),:);
    
    
    plot([1:nCohDisc]+xoff(ii), (grpTbl.mean_mean_accY(grpTbl.d5rt == ii)), '-', 'LineWidth', 2, 'Color', cols(cIdx(ii+1),:));
    
    allY = [allY; e.YData(:)];
    
end
set(gca, 'TickDir', 'out', 'LineWidth', 1)
xlim([.5, nCohDisc+.5])

yticks(ylim)
xticks([]);

yticks([ceil(min(allY)*100)/100, floor(max(allY)*100)/100])

title('Target \times RT')
xlabel('Target Coherence')

ylabel('relative accuracy')
yticks(ylim)
xticks([]);









% targ x RT (STATIC SIM)
grpTarg = grpstats(static_tbl, {'rt5', 'targ7'});
figM = nexttile; hold on;
for ii = 1:5
    plot((grpTarg.mean_acc(grpTarg.rt5 == ii)),...
        '-o', 'LineWidth', 2, 'Color', cols(cIdx(ii+1),:), 'MarkerFaceColor', cols(cIdx(ii+1),:));
end
set(gca, 'TickDir', 'out', 'LineWidth', 1);
title('static sim')
yticks(ylim)
xticks([]);
xlabel('target coherence')
ylabel('relative accuracy')
xlim([.5, nCohDisc+.5])




% allign axes
% lims = [min([figP.YAxis(1).Limits, figM.YAxis(1).Limits]), max([figP.YAxis(1).Limits, figM.YAxis(1).Limits])]
lims = [min([figP.YAxis(1).Limits, figM.YAxis(1).Limits]), max([figP.YAxis(1).Limits, figM.YAxis(1).Limits])];


figP.YAxis.Limits = lims; figP.YAxis.TickValues = lims;
figM.YAxis.Limits = lims; figM.YAxis.TickValues = lims; 











% ============== targ dynamics RT ==============
grpTbl = grpstats(rt_tbl, {'pt', 'acc', 'targ'}, {'mean', @skewness, @kurtosis});
grpTbl = grpstats(grpTbl, {'acc', 'targ'}, {'mean','sem'});




rt_cIdx = [180, 100];

xoff = [-.15, .15];

figP=nexttile([1,1]); hold on;
allY = [];
for ii = 0:1
    
    
        e=errorbar([1:nCohDisc]  + xoff(ii+1), (exp(grpTbl.mean_mean_rt(grpTbl.acc == ii)) + mean(all_ndt)), exp(grpTbl.mean_mean_rt(grpTbl.acc == ii)).*grpTbl.sem_mean_crt(grpTbl.acc == ii),...
            'o', 'LineWidth', 1);
    
%     e=errorbar([1:nCohDisc]  + xoff(ii+1), grpTbl.mean_skewness_rawRt(grpTbl.acc == ii), grpTbl.sem_skewness_rawRt(grpTbl.acc == ii),...
%         'o', 'LineWidth', 1);
    
    e.Color = cols(rt_cIdx(ii+1),:);
    e.MarkerFaceColor = cols(rt_cIdx(ii+1),:);
    
    
    plot([1:nCohDisc] + xoff(ii+1), (exp(grpTbl.mean_mean_rtY(grpTbl.acc == ii))+ mean(all_ndt)), '-', 'LineWidth', 2, 'Color', cols(rt_cIdx(ii+1),:));

    allY = [allY; e.YData(:)];
    
    
end
set(gca, 'TickDir', 'out', 'LineWidth', 1)
xlim([.5, nCohDisc+.5])
xticks([]);
title('Target \times Accuracy')

yticks([ceil(min(allY)*100)/100, floor(max(allY)*100)/100])

ylabel('RT')
yticks(ylim)
xticks([]);







% targ x Acc (STATIC SIM)
grpTarg = grpstats(static_tbl, {'acc', 'targ7'});
figM=nexttile; hold on;
for ii = 0:1
    plot((grpTarg.mean_rt(grpTarg.acc == ii)),...
        '-o', 'LineWidth', 2, 'Color', cols(rt_cIdx(ii+1),:), 'MarkerFaceColor', cols(rt_cIdx(ii+1),:));
end
set(gca, 'TickDir', 'out', 'LineWidth', 1);
title('static sim')
yticks(ylim)
xticks([]);
xlabel('target coherence')
ylabel('RT')
xlim([.5, nCohDisc+.5])


% allign axes
lims = [min([figP.YAxis(1).Limits, figM.YAxis(1).Limits]), max([figP.YAxis(1).Limits, figM.YAxis(1).Limits])];

figP.YAxis.Limits = lims; figP.YAxis.TickValues = lims;
figM.YAxis.Limits = lims; figM.YAxis.TickValues = lims; 





% ============ Distractor









% ============== dist dynamics acc ==============
grpTbl = grpstats(acc_tbl, {'pt', 'd5rt', 'dist'});
grpTbl = grpstats(grpTbl, {'d5rt', 'dist'}, {'mean','sem'});

figP=nexttile([1,1]); hold on;

xoff = linspace(-.2,.2,5);

allY = [];
for ii = 1:5
    
    e=errorbar([1:nCohDisc]+xoff(ii), center(grpTbl.mean_mean_acc(grpTbl.d5rt == ii)), grpTbl.sem_mean_cacc(grpTbl.d5rt == ii),...
        'o', 'LineWidth', 1);
    e.Color = cols(cIdx(ii+1),:);
    e.MarkerFaceColor = cols(cIdx(ii+1),:);
    
    plot([1:nCohDisc]+xoff(ii), center(grpTbl.mean_mean_accY(grpTbl.d5rt == ii)), '-', 'LineWidth', 2, 'Color', cols(cIdx(ii+1),:));
    
    allY = [allY; e.YData(:)];
    
end
set(gca, 'TickDir', 'out', 'LineWidth', 1)
xlim([.5, nCohDisc+.5])
xticks([]);

% legend('Fastest RT','','','','Slowest RT', 'Location', 'southeast')
ylabel('Relative Accuracy')
title('Distractor \times RT')
xlabel('distractor congruence')

yticks(ylim)
xticks([]);
yline(0, '-k', 'LineWidth', 1);




% dist X RT  (STATIC SIM)
grpDist = grpstats(static_tbl, {'rt5', 'dist7'});
figM=nexttile; hold on;
for ii = 1:5
    plot(center(grpDist.mean_acc(grpDist.rt5 == ii)),...
        '-o', 'LineWidth', 2, 'Color', cols(cIdx(ii+1),:), 'MarkerFaceColor', cols(cIdx(ii+1),:));
end

set(gca, 'TickDir', 'out', 'LineWidth', 1);
title('static sim')
yticks(ylim)
xticks([]);
xlabel('distractor congruence')
ylabel('relative accuracy')
xlim([.5, nCohDisc+.5])

yline(0, '-k', 'LineWidth', 1);


% allign axes
lims = [min([figP.YAxis(1).Limits, figM.YAxis(1).Limits]), max([figP.YAxis(1).Limits, figM.YAxis(1).Limits])];

figP.YAxis.Limits = lims; figP.YAxis.TickValues = lims;
figM.YAxis.Limits = lims; figM.YAxis.TickValues = lims; 








% ============== dist dynamics RT ==============
grpTbl = grpstats(rt_tbl, {'pt', 'acc', 'dist'});
grpTbl = grpstats(grpTbl, {'acc', 'dist'}, {'mean','sem'});

xoff = [-.15, .15];


figP=nexttile([1,1]); hold on;
allY = [];
for ii = 0:1
    
    
    e=errorbar([1:nCohDisc] + xoff(ii+1), (exp(grpTbl.mean_mean_rt(grpTbl.acc == ii)) + mean(all_ndt)), exp(grpTbl.mean_mean_rt(grpTbl.acc == ii)).*grpTbl.sem_mean_crt(grpTbl.acc == ii),...
        'o', 'LineWidth', 1);
    e.Color = cols(rt_cIdx(ii+1),:);
    e.MarkerFaceColor = cols(rt_cIdx(ii+1),:);
    
    
    plot([1:nCohDisc] + xoff(ii+1), (exp(grpTbl.mean_mean_rtY(grpTbl.acc == ii))+ mean(all_ndt)), '-', 'LineWidth', 2, 'Color', cols(rt_cIdx(ii+1),:));
    
    allY = [allY; e.YData(:)];
    
end
set(gca, 'TickDir', 'out', 'LineWidth', 1)
xlim([.5, nCohDisc+.5])
xticks([]);
xlabel('Distractor Congruence')

yticks([ceil(min(allY)*100)/100, floor(max(allY)*100)/100])

ylabel('RT')
title('Distractor \times Accuracy')
yticks(ylim)
xticks([]);






% dist x accuracy (STATIC SIM)
grpDist = grpstats(static_tbl, {'acc', 'dist7'});
figM=nexttile; hold on;
for ii = 0:1
    plot((grpDist.mean_rt(grpDist.acc == ii)),...
        '-o', 'LineWidth', 2, 'Color', cols(rt_cIdx(ii+1),:), 'MarkerFaceColor', cols(rt_cIdx(ii+1),:));
end

set(gca, 'TickDir', 'out', 'LineWidth', 1);
title('static sim')
yticks(ylim)
xticks([]);
xlabel('distractor congruence')
ylabel('RT')
xlim([.5, nCohDisc+.5])




% allign axes
lims = [min([figP.YAxis(1).Limits, figM.YAxis(1).Limits]), max([figP.YAxis(1).Limits, figM.YAxis(1).Limits])];

figP.YAxis.Limits = lims; figP.YAxis.TickValues = lims;
figM.YAxis.Limits = lims; figM.YAxis.TickValues = lims; 





%% ~~~~ Updated FFIc dynamic ~~~~~~~

nReals = 10000;          % how many simulated trials?



% ===== simulation parameters  =====
threshType = 'race';    % which threshold {race, diff}
doRelu = 1;             % force positive accumulators?

nChoices = 2;           % how many accumulators?
nDt    = 2000;          % max number of timesteps
dt = .001;              % simulation dt

% ===== task parameters =====

nTarg = 11;
nDist = 11;
targRange = linspace(.6,.95,nTarg);
distRange = linspace(-1,1,nDist);

[targ, dist] = meshgrid(targRange, distRange);
targ = tanh(1.5*targ(:))./tanh(1.5);
ntarg = tanh(1.5*(1-targ(:)))./tanh(1.5);
dist = tanh(1.35*dist(:))./tanh(1.35);

nConds = length(targ);





% ========= run COLOR dynamic sim




% decision dynamics

decay = 0;          % accumulator-specific decay
cross = 0;          % accumulator competition
ffinh = 1;          % feedforward inhibition
sd    = 1;          % diffusion noise
lapse = 0;          % lapse rate


A = cross*ones(nChoices) + (decay - cross)*eye(nChoices)
B = [...
    1, -ffinh;...
    -ffinh, 1;...
    ]


% gain dynamics
targ0       = 1;
targAim     = 6;
targGain    = 3;
 

dist0       = 1.25;
distAim     = 0;
distGain    = 4.5;

salSD       = 0.1;

 % non-decision time
t0 = .300;         

% decision threshold
thresh = max(.01, linspace(1.5, .1, nDt))';         




% a =  ones(2)*-25
% B = [...
%     1, 0;...
%     0, 1;...
%     ]




sddt = sd*sqrt(dt);
salSDdt = salSD*sqrt(dt);


% run sims  ======================

clear meanRT Accuracy static_allRT static_allChoice

tic
parfor cc = 1:nConds
    
    x = zeros(nChoices,nReals);
    all_x = nan(nDt,nChoices,nReals);
    
    %     targ0       = 2;
    %     targAim     = 5;
    %     targGain    = 1.5;
    %
    %     dist0       = 1;
    %     distAim     = 0;
    %     distGain    = 1.5;
    
    targSal = targ0;
    distSal = dist0;
    

    I = [...
        targ(cc),  abs(dist(cc))*(dist(cc)>0);...
        ntarg(cc), abs(dist(cc))*(dist(cc)<0);...
        ];
    
    
    for tt = 1:nDt
        
        % update gain
        targSal = targSal + targGain*(targAim - targSal)*dt + randn*salSDdt;
        distSal = distSal + distGain*(distAim - distSal)*dt + randn*salSDdt;
        
        % update model
        all_x(tt,:,:) = x;
        x = x + A*(x*dt) + B*(I*[targSal;distSal]*dt + randn(size(x))*sddt);
        
        % threshold
        if doRelu
            x = max(x,0);
        end
        
        % break if all accumulators reach threshold
        switch threshType
            
            case 'race'
                
                if ~any(max(x)<thresh(tt))
                    break;
                end
                
            case 'diff'
                
                if ~any(abs(diff(x))<thresh(tt))
                    break;
                end
                
        end
        
    end
    
    [RT, Choice] = deal(nan(nReals, 1));
    
    for rr = 1:nReals
        try
            switch threshType
                case 'race'
                    RT(rr)      = find(max(all_x(:,:,rr),[],2) >= thresh, 1);
                case 'diff'
                    RT(rr)      = find(abs(diff(all_x(:,:,rr),[],2)) >= thresh, 1);
            end
            
            [~,Choice(rr)]  = max(all_x(RT(rr),:,rr));
            
        catch
            RT(rr) = nan;
            Choice(rr) = 0;
        end
    end
    
    stdRT(cc) = nanstd(RT(Choice == 1)*dt + t0);
    meanRT(cc) = nanmean(RT(Choice == 1)*dt + t0);
    Accuracy(cc) = (1-lapse)*nanmean(Choice == 1) + lapse*.5;
    
    static_allRT{cc} = RT;
    static_allChoice{cc} = Choice;
    
end
toc


% get mean
rStdRTStatic = reshape(stdRT,[nTarg, nDist]);
rRTStatic = reshape(meanRT,[nTarg, nDist]);
rAccStatic = reshape(Accuracy,[nTarg, nDist]);




% make table
[t_dist, t_targ, t_acc, t_rt] = deal([]);
for cc = 1:nConds
    
    t_acc  = [t_acc; static_allChoice{cc}==1];
    t_rt   = [t_rt; static_allRT{cc}*dt + t0];
    t_targ = [t_targ; ones(nReals,1)*targ(cc)];
    t_dist = [t_dist; ones(nReals,1)*dist(cc)];
    
end



static_tbl = table;
static_tbl.acc = t_acc;
static_tbl.rt = t_rt;
static_tbl.lrt = log(t_rt);
static_tbl.crt = t_rt - nanmean(t_rt);
static_tbl.rt5 = discretize(tiedrank(t_rt),5);
static_tbl.targ = t_targ;
static_tbl.targ7 = discretize(tiedrank(t_targ),7);
static_tbl.dist = t_dist;
static_tbl.dist7 = discretize(tiedrank(t_dist),7);



% regression


rt_mdl = fitlm(static_tbl, 'lrt ~ targ*dist', 'Exclude', static_tbl.acc == 0)
acc_mdl = fitglm(static_tbl, 'acc ~ targ*dist', 'Distribution', 'binomial')
dyn_mdl = fitglm(static_tbl, 'acc ~ (targ + dist)*crt', 'Distribution', 'binomial')











% ========= run MOTION simulation




% decision dynamics

decay = 0;          % accumulator-specific decay
cross = 0;          % accumulator competition
ffinh = 1;          % feedforward inhibition
sd    = 1;          % diffusion noise
lapse = 0;          % lapse rate


A = cross*ones(nChoices) + (decay - cross)*eye(nChoices)
B = [...
    1, -ffinh;...
    -ffinh, 1;...
    ]


% gain dynamics
targ0       = 5;
targAim     = 5;
targGain    = 4;
 

dist0       = 0;
distAim     = 0;
distGain    = 4;

salSD       = 0.1;

 % non-decision time
t0 = .200;         

% decision threshold
thresh = max(.01, linspace(1.5, .1, nDt))';         




% a =  ones(2)*-25
% B = [...
%     1, 0;...
%     0, 1;...
%     ]




sddt = sd*sqrt(dt);
salSDdt = salSD*sqrt(dt);


% run sims  ======================

clear meanRT Accuracy static_allRT static_allChoice

tic
parfor cc = 1:nConds
    
    x = zeros(nChoices,nReals);
    all_x = nan(nDt,nChoices,nReals);
    
    %     targ0       = 2;
    %     targAim     = 5;
    %     targGain    = 1.5;
    %
    %     dist0       = 1;
    %     distAim     = 0;
    %     distGain    = 1.5;
    
    targSal = targ0;
    distSal = dist0;
    

    I = [...
        targ(cc),  abs(dist(cc))*(dist(cc)>0);...
        ntarg(cc), abs(dist(cc))*(dist(cc)<0);...
        ];
    
    
    for tt = 1:nDt
        
        % update gain
        targSal = targSal + targGain*(targAim - targSal)*dt + randn*salSDdt;
        distSal = distSal + distGain*(distAim - distSal)*dt + randn*salSDdt;
        
        % update model
        all_x(tt,:,:) = x;
        x = x + A*(x*dt) + B*(I*[targSal;distSal]*dt + randn(size(x))*sddt);
        
        % threshold
        if doRelu
            x = max(x,0);
        end
        
        % break if all accumulators reach threshold
        switch threshType
            
            case 'race'
                
                if ~any(max(x)<thresh(tt))
                    break;
                end
                
            case 'diff'
                
                if ~any(abs(diff(x))<thresh(tt))
                    break;
                end
                
        end
        
    end
    
    [RT, Choice] = deal(nan(nReals, 1));
    
    for rr = 1:nReals
        try
            switch threshType
                case 'race'
                    RT(rr)      = find(max(all_x(:,:,rr),[],2) >= thresh, 1);
                case 'diff'
                    RT(rr)      = find(abs(diff(all_x(:,:,rr),[],2)) >= thresh, 1);
            end
            
            [~,Choice(rr)]  = max(all_x(RT(rr),:,rr));
            
        catch
            RT(rr) = nan;
            Choice(rr) = 0;
        end
    end
    
    stdRT(cc) = nanstd(RT(Choice == 1)*dt + t0);
    meanRT(cc) = nanmean(RT(Choice == 1)*dt + t0);
    Accuracy(cc) = (1-lapse)*nanmean(Choice == 1) + lapse*.5;
    
    static_allRT{cc} = RT;
    static_allChoice{cc} = Choice;
    
end
toc



% get mean
rRTMot = reshape(meanRT,[nTarg, nDist]);
rAccMot = reshape(Accuracy,[nTarg, nDist]);







 % PLOT within-trial  



 

% ===== start plot
figure('units','inch','position',[0,0,1.5*7.5,1.5*7.5]); hold on;
tiledlayout(4,4, 'TileSpacing', 'compact', 'Padding', 'compact');


 
 
 
 
 
 


load('batlow.mat');
cols = batlow;
cIdx = round(linspace(1, length(cols), 7));

rtacc_cols = {[36, 123, 160]/255, [242, 95, 92]/255};




mdl = wn_col;


[rt, crt, rt_pred, rt_targ, rt_dist, rt_res, rt_pt,...
    acc, cacc, acc_pred, acc_targ, acc_dist, acc_res, acc_pt,...
    all_ndt,...
    ] = deal([]);

pt_count = 1;

for mm = 1:length(col)
    
    
    data = mdl{mm}.data;
    results = mdl{mm}.results;
    
    param = results.E';
    ndt = ...
        results.link{ismember(results.paramName, 'rt_{ndt}')} ...
        (param(:, ismember(results.paramName, 'rt_{ndt}')));
    
    all_ndt = [all_ndt; ndt];
    
    for tt = 1:size(data,2)
        
        %                 if mm == 1 && tt == size(data,2)
        %                     continue
        %                 end
        
        % simulation parameters
        [~, out]=eval([results.likfun, '(param(tt,:), data(tt), [], [], 0, [])']);
        
        
        prt         = log(out.rt.rt(out.rt.rt_sel) - ndt(tt));
        prt(out.rt.rt(out.rt.rt_sel) <= ndt(tt)) = nan;
        rt          = [rt; prt];
        crt          = [crt; prt - nanmean(prt)];
        
        
        rt_pred     = [rt_pred; out.rt.rt_yhat];
        
        rt_targ     = [rt_targ; out.data.d11Targ0(out.rt.rt_sel)];
        rt_dist     = [rt_dist; 12-out.data.d11Dist0(out.rt.rt_sel)];
        %         rt_targ    = [rt_targ; out.data.d11Targ0(out.rt.rt_sel).*out.data.sgnResp(out.rt.rt_sel)];
        %         rt_dist    = [rt_dist; (12-out.data.d11Dist0(out.rt.rt_sel)).*out.data.sgnResp(out.rt.rt_sel)];
        
        rt_res      = [rt_res; prt - out.rt.rt_yhat];
        
        rt_pt       = [rt_pt; out.rt.rt_X(:,1).*pt_count];
        
        
        %         acc         = [acc; out.acc.acc(out.acc.acc_sel)];
        %         acc_pred    = [acc_pred; out.acc.acc_yhat];
        %         acc_targ    = [acc_targ; out.data.d11Targ0(out.acc.acc_sel).*out.data.sgnResp(out.acc.acc_sel)];
        %         acc_dist    = [acc_dist; (12-out.data.d11Dist0(out.acc.acc_sel)).*out.data.sgnResp(out.acc.acc_sel)];
        %         acc_pt      = [acc_pt; out.acc.acc_X(:,1).*pt_count];
        
        acc         = [acc; out.data.acc0(out.acc.acc_sel)];
        cacc        = [cacc; center(out.data.acc0(out.acc.acc_sel))];
        
        yhat        = [out.acc.acc_yhat];
        yhat(out.data.corrResp(out.acc.acc_sel)==0) = 1-yhat(out.data.corrResp(out.acc.acc_sel)==0);
        acc_pred    = [acc_pred; yhat];
        
        
        acc_targ    = [acc_targ; out.data.d11Targ0(out.acc.acc_sel)];
        acc_dist    = [acc_dist; (12-out.data.d11Dist0(out.acc.acc_sel))];
        acc_pt      = [acc_pt; out.acc.acc_X(:,1).*pt_count];
        
        pt_count = pt_count+1;
        
    end
    
    
end


cols = {[36, 123, 160]/255, [242, 95, 92]/255};



rt_tbl = array2table([rt, crt, rt_pred, rt_targ, rt_dist, rt_res, rt_pt ], 'VariableNames', ...
    {'rt', 'crt', 'rtY', 'targAll', 'distAll', 'res', 'pt'});

acc_tbl = array2table([acc, cacc, acc_pred, acc_targ, acc_dist, acc_pt], 'VariableNames', ...
    {'acc', 'cacc', 'accY', 'targAll', 'distAll', 'pt'});










% ===== MOTION


mdl = wn_mot;
% mdl = exp23;


[rt, crt, d5rt,rawRt, rt_pred, rt_targ, rt_dist, rt_acc, rt_pt,...
    acc, cacc, acc_pred, acc_targ, acc_dist, acc_rt, acc_pt,...
    all_ndtMot, rt_targ_all, rt_dist_all, acc_targ_all, acc_dist_all, ...
    ] = deal([]);

prevCor = [];

pt_count = 1;

nCohDisc = 7;


for mm = 1:3
    
    
    data = mdl{mm}.data;
    results = mdl{mm}.results;
    
    param = results.E';
    ndt = ...
        results.link{ismember(results.paramName, 'rt_{ndt}')} ...
        (param(:, ismember(results.paramName, 'rt_{ndt}')));
    
    all_ndtMot = [all_ndtMot; ndt];
    
    for tt = 1:size(data,2)
        
        %                 if mm == 1 && tt == size(data,2)
        %                     continue
        %                 end
        
        % simulation parameters
        [~, out]=eval([results.likfun, '(transpose(param(tt,:)), data(tt), [], [], 0, [])']);
        
        
        % rt
        prt         = log(out.rt.rt(out.rt.rt_sel) - ndt(tt));
        prt(out.rt.rt(out.rt.rt_sel) <= ndt(tt)) = nan;
        rt          = [rt; prt];
        crt          = [crt; prt - nanmean(prt)];
        rawRt       = [rawRt; out.rt.rt(out.rt.rt_sel)];
        
        
        rt_pred     = [rt_pred; out.rt.rt_yhat];
        
        rt_targ     = [acc_targ; discretize(tiedrank(out.data.targ0(out.rt.rt_sel)),nCohDisc)];
        rt_dist     = [acc_dist; discretize(tiedrank(out.data.dist0(out.rt.rt_sel)),nCohDisc)];
        
        
        rt_targ_all     = [rt_targ_all; out.data.d11Targ0(out.rt.rt_sel)];
        rt_dist_all     = [rt_dist_all; 12-out.data.d11Dist0(out.rt.rt_sel)];
        
        
        rt_acc      = [rt_acc; out.data.acc0(out.rt.rt_sel)];
        
        rt_pt       = [rt_pt; out.rt.rt_X(:,1).*pt_count];
        
        
        % acc
        d5rt          = [d5rt; discretize(tiedrank(out.data.rt0(out.acc.acc_sel)),5)];
        
        acc         = [acc; out.data.acc0(out.acc.acc_sel)];
        cacc        = [cacc; center(out.data.acc0(out.acc.acc_sel))];
        
        yhat        = [out.acc.acc_yhat];
        yhat(out.data.corrResp(out.acc.acc_sel)==0) = 1-yhat(out.data.corrResp(out.acc.acc_sel)==0);
        acc_pred    = [acc_pred; yhat];
        
        
        acc_targ     = [acc_targ; discretize(tiedrank(out.data.targ0(out.acc.acc_sel)),nCohDisc)];
        acc_dist     = [acc_dist; discretize(tiedrank(out.data.dist0(out.acc.acc_sel)),nCohDisc)];
        
        acc_targ_all     = [acc_targ_all; out.data.d11Targ0(out.acc.acc_sel)];
        acc_dist_all     = [acc_dist_all; 12-out.data.d11Dist0(out.acc.acc_sel)];
        
        acc_rt          = [acc_rt; out.data.rt0(out.acc.acc_sel)];
        acc_pt      = [acc_pt; out.acc.acc_X(:,1).*pt_count];
        
        pt_count = pt_count+1;
        
    end
    
    
end



rtMot_tbl = array2table([rt, crt, rt_pred, rt_targ_all, rt_dist_all, rt_acc, rt_pt ], 'VariableNames', ...
    {'rt', 'crt', 'rtY',  'targAll', 'distAll', 'acc', 'pt'});

accMot_tbl = array2table([acc, cacc, acc_pred, acc_rt, d5rt, acc_targ_all, acc_dist_all, acc_pt], 'VariableNames', ...
    {'acc', 'cacc', 'accY', 'rt', 'd5rt',  'targAll', 'distAll', 'pt'});










% acc_tbl = array2table([acc, acc_pred, rt, d5rt, acc_targ, acc_dist, acc_pt], 'VariableNames', ...
%     {'acc', 'accY','rt', 'd5rt', 'targ', 'dist' 'pt'});


npt = [size(col{1}.data,2), size(col{2}.data,2), size(col{3}.data,2)];
nptSel = [ones(npt(1),1)*1; ones(npt(2),1)*2; ones(npt(3),1)*3];











% =================== TARGET


figP = nexttile([1,1]); hold on;

% ============== targ acc ==============
grpTbl = grpstats(acc_tbl, {'pt','targAll'});
grpTbl = grpstats(grpTbl, {'targAll'}, {'mean','sem'});

yyaxis('left');

% plot([1:length(grpTbl.mean_mean_acc)],(grpTbl.mean_mean_acc),'ob', 'LineWidth', .5);
errorbar(1.15:11.15, grpTbl.mean_mean_acc, grpTbl.sem_mean_cacc,...
    'o', 'color', rtacc_cols{1}, 'MarkerFaceColor', rtacc_cols{1}, 'LineWidth', 1);

plot(1.15:11.15, grpTbl.mean_mean_accY, '-', 'color', rtacc_cols{1}, 'LineWidth', 2);
set(gca, 'TickDir', 'out', 'LineWidth', 1)

xticks([]);

% ============== targ rt ==============
grpTbl = grpstats(rt_tbl, {'pt','targAll'}, {'mean'});
grpTbl = grpstats(grpTbl, {'targAll'}, {'mean','sem'});

yyaxis('right');

% plot(exp(grpTbl.mean_mean_rt) + mean(all_ndt), 'ob', 'LineWidth', .5);
figP_RT = errorbar(.85:10.85, exp(grpTbl.mean_mean_rt) + mean(all_ndt), exp(grpTbl.mean_mean_rt).*(grpTbl.sem_mean_crt),...
    'o', 'color', rtacc_cols{2}, 'MarkerFaceColor', rtacc_cols{2}, 'LineWidth', 1);

plot(.85:10.85, exp(grpTbl.mean_mean_rtY) + mean(all_ndt), '-', 'color', rtacc_cols{2}, 'LineWidth', 2);

set(gca, 'TickDir', 'out', 'LineWidth', 1)

yticks(ylim)

title('Target (Attend-Color)')

xlim([.5 11.5])
xticks([])
xticklabels({'low coh', 'high coh'})

xlabel('target coherence')


colTargRTY = ylim;
yticks(colTargRTY)
xticks([]);





% plot target DDM color =========


% target
figM = nexttile([1,1]); hold on;

yyaxis('left');
plot(.85:10.85, mean(rAccStatic,1), '-o', 'color', rtacc_cols{1}, 'MarkerFaceColor', rtacc_cols{1}, 'LineWidth', 2);
ylabel('accuracy')

xticks([]);



yyaxis('right');
plot(.85:10.85, mean(rRTStatic,1), '-o', 'color', rtacc_cols{2}, 'MarkerFaceColor', rtacc_cols{2}, 'LineWidth', 2);
ylabel('RT')

set(gca, 'LineWidth', 1, 'TickDir', 'out');
title('DDM sim')
xlabel('target coherence')

xticks([]);
xlim([.5 11.5])


% allign axes
accLim = [min([figP.YAxis(1).Limits, figM.YAxis(1).Limits])-.01, max([figP.YAxis(1).Limits, figM.YAxis(1).Limits])+.01]
rtLim = [min([figP.YAxis(2).Limits, figM.YAxis(2).Limits]), max([figP.YAxis(2).Limits, figM.YAxis(2).Limits])]

figP.YAxis(1).Limits = accLim; figP.YAxis(1).TickValues = accLim;
figM.YAxis(1).Limits = accLim; figM.YAxis(1).TickValues = accLim;

figP.YAxis(2).Limits = rtLim; figP.YAxis(2).TickValues = rtLim; 
figM.YAxis(2).Limits = rtLim; figM.YAxis(2).TickValues = rtLim;





% Attend-MOTION (target)




figP = nexttile([1,1]); hold on;

% ============== targ acc ==============
grpTbl = grpstats(accMot_tbl, {'pt','targAll'});
grpTbl = grpstats(grpTbl, {'targAll'}, {'mean','sem'});

yyaxis('left');

% plot([1:length(grpTbl.mean_mean_acc)],(grpTbl.mean_mean_acc),'ob', 'LineWidth', .5);
ep = errorbar(1.15:11.15, grpTbl.mean_mean_acc, grpTbl.sem_mean_cacc,...
    'o', 'color', rtacc_cols{1}, 'MarkerFaceColor', rtacc_cols{1}, 'LineWidth', 1);

plot(1.15:11.15, grpTbl.mean_mean_accY, '-', 'color', rtacc_cols{1}, 'LineWidth', 2);
set(gca, 'TickDir', 'out', 'LineWidth', 1)
yticks([floor(min(ep.YData)*100)/100, ceil(max(ep.YData)*100)/100])


motTargAccY = ylim;
yticks(motTargAccY)
xticks([]);




% ============== targ rt ==============
grpTbl = grpstats(rtMot_tbl, {'pt','targAll'}, {'mean'});
grpTbl = grpstats(grpTbl, {'targAll'}, {'mean','sem'});

yyaxis('right');

% plot(exp(grpTbl.mean_mean_rt) + mean(all_ndt), 'ob', 'LineWidth', .5);
ep=errorbar(.85:10.85, exp(grpTbl.mean_mean_rt) + mean(all_ndtMot), exp(grpTbl.mean_mean_rt).*(grpTbl.sem_mean_crt),...
    'o', 'color', rtacc_cols{2}, 'MarkerFaceColor', rtacc_cols{2}, 'LineWidth', 1);

plot(.85:10.85, exp(grpTbl.mean_mean_rtY) + mean(all_ndtMot), '-', 'color', rtacc_cols{2}, 'LineWidth', 2);

set(gca, 'TickDir', 'out', 'LineWidth', 1)
yticks([floor(min(ep.YData)*100)/100, ceil(max(ep.YData)*100)/100])


title('Target (Attend-Motion)')

xlim([.5 11.5])
xticks([])
xticklabels({'low coh', 'high coh'})

xlabel('target coherence')


motTargRTY = ylim;
yticks(motTargRTY)
xticks([]);




% plot target DDM motion =========


% target
figM = nexttile([1,1]); hold on;

yyaxis('left');
plot(.85:10.85, mean(rAccMot,1), '-o', 'color', rtacc_cols{1}, 'MarkerFaceColor', rtacc_cols{1}, 'LineWidth', 2);
ylabel('accuracy')
xticks([]);
yticks([]);

xlim([.5 11.5])


% ylim(motTargAccY);
yticks(ylim)
xticks([]);






yyaxis('right');
plot(.85:10.85, mean(rRTMot,1), '-o', 'color', rtacc_cols{2}, 'MarkerFaceColor', rtacc_cols{2}, 'LineWidth', 2);
ylabel('RT')

xlim([.5 11.5])
set(gca, 'LineWidth', 1, 'TickDir', 'out');
title('DDM sim')
xlabel('target coherence')

% ylim(motTargRTY);
yticks(ylim)
xticks([]);



% allign axes
accLim = [min([figP.YAxis(1).Limits, figM.YAxis(1).Limits]), max([figP.YAxis(1).Limits, figM.YAxis(1).Limits])]
rtLim = [min([figP.YAxis(2).Limits, figM.YAxis(2).Limits]), max([figP.YAxis(2).Limits, figM.YAxis(2).Limits])]

figP.YAxis(1).Limits = accLim; figP.YAxis(1).TickValues = accLim;
figM.YAxis(1).Limits = accLim; figM.YAxis(1).TickValues = accLim;

figP.YAxis(2).Limits = rtLim; figP.YAxis(2).TickValues = rtLim; 
figM.YAxis(2).Limits = rtLim; figM.YAxis(2).TickValues = rtLim;











% distractor


figP = nexttile([1,1]); hold on;

% ============== dist acc ==============
yyaxis('left');

grpTbl = grpstats(acc_tbl, {'pt','distAll'});
grpTbl = grpstats(grpTbl, {'distAll'}, {'mean','sem'});

% plot([1:length(grpTbl.mean_mean_acc)], (grpTbl.mean_mean_acc),'ok', 'LineWidth', .5);
ep=errorbar(1.15:11.15, grpTbl.mean_mean_acc, grpTbl.sem_mean_cacc,...
    'o', 'color', rtacc_cols{1}, 'MarkerFaceColor', rtacc_cols{1}, 'LineWidth', 1);
plot(1.15:11.15, grpTbl.mean_mean_accY, '-', 'color', rtacc_cols{1}, 'LineWidth', 2);
set(gca, 'TickDir', 'out', 'LineWidth', 1)
yticks([floor(min(ep.YData)*100)/100, ceil(max(ep.YData)*100)/100])
xticks([])


motDistAccY = ylim;
yticks(motDistAccY)
xticks([]);


% ============== dist rt ==============
yyaxis('right');

grpTbl = grpstats(rt_tbl, {'pt','distAll'}, {'mean'});
grpTbl = grpstats(grpTbl, {'distAll'}, {'mean','sem'});



% plot(exp(grpTbl.mean_mean_rt) + mean(all_ndt), 'ok', 'LineWidth', .5);
ep = errorbar(.85:10.85, exp(grpTbl.mean_mean_rt) + mean(all_ndt), exp(grpTbl.mean_mean_rt).*(grpTbl.sem_mean_crt),...
    'o', 'color', rtacc_cols{2}, 'MarkerFaceColor', rtacc_cols{2}, 'LineWidth', 1);
plot(.85:10.85, exp(grpTbl.mean_mean_rtY)+ mean(all_ndt), '-', 'color', rtacc_cols{2}, 'LineWidth', 2);
set(gca, 'TickDir', 'out', 'LineWidth', 1)
title('Distractor (Attend-Color)')
% ylim(targRange);

yticks([floor(min(ep.YData)*100)/100, ceil(max(ep.YData)*100)/100])
xlim([.5 11.5])
xticks([])
xlabel('distractor congruence')


motDistRTY = ylim;
yticks(ylim)
xticks([]);












% plot target DDM color =========


figM = nexttile([1,1]); hold on;

yyaxis('left');
plot(.85:10.85, mean(rAccStatic,2), '-o', 'color', rtacc_cols{1}, 'MarkerFaceColor', rtacc_cols{1}, 'LineWidth', 2);
ylabel('accuracy')


% ylim(motDistAccY);
yticks(ylim)
xticks([]);


xlim([.5 11.5])

yyaxis('right');
plot(.85:10.85, mean(rRTStatic,2), '-o', 'color', rtacc_cols{2}, 'MarkerFaceColor', rtacc_cols{2}, 'LineWidth', 2);
ylabel('RT')

set(gca, 'LineWidth', 1, 'TickDir', 'out');
title('DDM sim')
xlabel('distractor congruence')


% ylim(motDistRTY);
yticks(ylim)
xticks([]);



% allign axes
accLim = [min([figP.YAxis(1).Limits, figM.YAxis(1).Limits]), max([figP.YAxis(1).Limits, figM.YAxis(1).Limits])]
rtLim = [min([figP.YAxis(2).Limits, figM.YAxis(2).Limits]), max([figP.YAxis(2).Limits, figM.YAxis(2).Limits])]

figP.YAxis(1).Limits = accLim; figP.YAxis(1).TickValues = accLim;
figM.YAxis(1).Limits = accLim; figM.YAxis(1).TickValues = accLim;

figP.YAxis(2).Limits = rtLim; figP.YAxis(2).TickValues = rtLim; 
figM.YAxis(2).Limits = rtLim; figM.YAxis(2).TickValues = rtLim;









% motion distractor


figP = nexttile([1,1]); hold on;

% ============== dist acc ==============
yyaxis('left');

grpTbl = grpstats(accMot_tbl, {'pt','distAll'});
grpTbl = grpstats(grpTbl, {'distAll'}, {'mean','sem'});

% plot([1:length(grpTbl.mean_mean_acc)], (grpTbl.mean_mean_acc),'ok', 'LineWidth', .5);
ep=errorbar(1.15:11.15, grpTbl.mean_mean_acc, grpTbl.sem_mean_cacc,...
    'o', 'color', rtacc_cols{1}, 'MarkerFaceColor', rtacc_cols{1}, 'LineWidth', 1);
plot(1.15:11.15, grpTbl.mean_mean_accY, '-', 'color', rtacc_cols{1}, 'LineWidth', 2);
set(gca, 'TickDir', 'out', 'LineWidth', 1)


ylim([.9, 1])
yticks(ylim)
xticks([]);


% ============== dist rt ==============
yyaxis('right');

grpTbl = grpstats(rtMot_tbl, {'pt','distAll'}, {'mean'});
grpTbl = grpstats(grpTbl, {'distAll'}, {'mean','sem'});



% plot(exp(grpTbl.mean_mean_rt) + mean(all_ndt), 'ok', 'LineWidth', .5);
ep = errorbar(.85:10.85, exp(grpTbl.mean_mean_rt) + mean(all_ndtMot), exp(grpTbl.mean_mean_rt).*(grpTbl.sem_mean_crt),...
    'o', 'color', rtacc_cols{2}, 'MarkerFaceColor', rtacc_cols{2}, 'LineWidth', 1);
plot(.85:10.85, exp(grpTbl.mean_mean_rtY)+ mean(all_ndtMot), '-', 'color', rtacc_cols{2}, 'LineWidth', 2);
set(gca, 'TickDir', 'out', 'LineWidth', 1)
title('Distractor (Attend-Motion)')
% ylim(targRange);

yticks([floor(min(ep.YData)*100)/100, ceil(max(ep.YData)*100)/100])
xlim([.5 11.5])
xticks([])
xlabel('distractor congruence')

ylim([.46, .52])
yticks(ylim)
xticks([]);






% plot target DDM color =========


figM = nexttile([1,1]); hold on;

yyaxis('left');
plot(.85:10.85, mean(rAccMot,2), '-o', 'color', rtacc_cols{1}, 'MarkerFaceColor', rtacc_cols{1}, 'LineWidth', 2);
ylabel('accuracy')

xlim([.5 11.5])
ylim([.9, 1])
yticks(ylim)
xticks([]);


yyaxis('right');
plot(.85:10.85, mean(rRTMot,2), '-o', 'color', rtacc_cols{2}, 'MarkerFaceColor', rtacc_cols{2}, 'LineWidth', 2);
ylabel('RT')

ylim([mean(ylim)-.03, mean(ylim)+.03])


set(gca, 'LineWidth', 1, 'TickDir', 'out');
title('DDM sim')
xlabel('distractor congruence')

ylim([.46, .52])
yticks(ylim)
xticks([]);




% allign axes
accLim = [min([figP.YAxis(1).Limits, figM.YAxis(1).Limits]), max([figP.YAxis(1).Limits, figM.YAxis(1).Limits])]
rtLim = [min([figP.YAxis(2).Limits, figM.YAxis(2).Limits]), max([figP.YAxis(2).Limits, figM.YAxis(2).Limits])]

figP.YAxis(1).Limits = accLim; figP.YAxis(1).TickValues = accLim;
figM.YAxis(1).Limits = accLim; figM.YAxis(1).TickValues = accLim;

figP.YAxis(2).Limits = rtLim; figP.YAxis(2).TickValues = rtLim; 
figM.YAxis(2).Limits = rtLim; figM.YAxis(2).TickValues = rtLim;









% ======= DYNAMICS





load('batlow.mat');
cols = batlow;
cIdx = round(linspace(1, length(cols), 7));

rtacc_cols = {[36, 123, 160]/255, [242, 95, 92]/255};









mdl = col;
% mdl = exp23;


[rt, rawRt, crt, d5rt, rt_pred, rt_targ, rt_dist, rt_acc, rt_pt,...
    acc, cacc, acc_pred, acc_targ, acc_dist, acc_rt, acc_pt,...
    all_ndt, rt_targ_all, rt_dist_all, acc_targ_all, acc_dist_all, ...
    ] = deal([]);

prevCor = [];

pt_count = 1;

nCohDisc = 7;


for mm = 1:3
    
    
    data = mdl{mm}.data;
    results = mdl{mm}.results;
    
    param = results.E';
    ndt = ...
        results.link{ismember(results.paramName, 'rt_{ndt}')} ...
        (param(:, ismember(results.paramName, 'rt_{ndt}')));
    
    all_ndt = [all_ndt; ndt];
    
    for tt = 1:size(data,2)
        
        %                 if mm == 1 && tt == size(data,2)
        %                     continue
        %                 end
        
        % simulation parameters
        [~, out]=eval([results.likfun, '(transpose(param(tt,:)), data(tt), [], [], 0, [])']);
        
        
        % rt
        prt         = log(out.rt.rt(out.rt.rt_sel) - ndt(tt));
        prt(out.rt.rt(out.rt.rt_sel) <= ndt(tt)) = nan;
        rt          = [rt; prt];
        crt         = [crt; prt - nanmean(prt)];
        rawRt       = [rawRt; out.rt.rt(out.rt.rt_sel)];

        
        rt_pred     = [rt_pred; out.rt.rt_yhat];
        
        rt_targ     = [acc_targ; discretize(tiedrank(out.data.targ0(out.rt.rt_sel)),nCohDisc)];
        rt_dist     = [acc_dist; discretize(tiedrank(out.data.dist0(out.rt.rt_sel)),nCohDisc)];
        
        
        rt_targ_all     = [rt_targ_all; out.data.d11Targ0(out.rt.rt_sel)];
        rt_dist_all     = [rt_dist_all; 12-out.data.d11Dist0(out.rt.rt_sel)];
        
        
        rt_acc      = [rt_acc; out.data.acc0(out.rt.rt_sel)];
        
        rt_pt       = [rt_pt; out.rt.rt_X(:,1).*pt_count];
        
        
        % acc
        d5rt          = [d5rt; discretize(tiedrank(out.data.rt0(out.acc.acc_sel)),5)];
        
        acc         = [acc; out.data.acc0(out.acc.acc_sel)];
        cacc        = [cacc; center(out.data.acc0(out.acc.acc_sel))];
        
        yhat        = [out.acc.acc_yhat];
        yhat(out.data.corrResp(out.acc.acc_sel)==0) = 1-yhat(out.data.corrResp(out.acc.acc_sel)==0);
        acc_pred    = [acc_pred; yhat];
        
        
        acc_targ     = [acc_targ; discretize(tiedrank(out.data.targ0(out.acc.acc_sel)),nCohDisc)];
        acc_dist     = [acc_dist; discretize(tiedrank(out.data.dist0(out.acc.acc_sel)),nCohDisc)];
        
        acc_targ_all     = [acc_targ_all; out.data.d11Targ0(out.acc.acc_sel)];
        acc_dist_all     = [acc_dist_all; 12-out.data.d11Dist0(out.acc.acc_sel)];
        
        acc_rt          = [acc_rt; out.data.rt0(out.acc.acc_sel)];
        acc_pt      = [acc_pt; out.acc.acc_X(:,1).*pt_count];
        
        pt_count = pt_count+1;
        
    end
    
    
end



rt_tbl = array2table([rt, rawRt, crt, rt_pred, rt_targ, rt_dist, rt_targ_all, rt_dist_all, rt_acc, rt_pt ], 'VariableNames', ...
    {'rt', 'rawRt', 'crt', 'rtY', 'targ', 'dist', 'targAll', 'distAll', 'acc', 'pt'});

acc_tbl = array2table([acc, cacc, acc_pred, acc_rt, d5rt, acc_targ, acc_dist,  acc_targ_all, acc_dist_all, acc_pt], 'VariableNames', ...
    {'acc', 'cacc', 'accY', 'rt', 'd5rt', 'targ', 'dist', 'targAll', 'distAll', 'pt'});









% ============== plot targ dynamics acc ==============
grpTbl = grpstats(acc_tbl, {'pt', 'd5rt', 'targ'});
grpTbl = grpstats(grpTbl, {'d5rt', 'targ'}, {'mean','sem'});


figP = nexttile([1,1]); hold on;

xoff = linspace(-.2,.2,5);
allY = [];
for ii = 1:5
    
    mAcc = mean(grpTbl.mean_mean_acc(grpTbl.d5rt == ii));
    
    e=errorbar([1:nCohDisc]+xoff(ii), (grpTbl.mean_mean_acc(grpTbl.d5rt == ii)), grpTbl.sem_mean_cacc(grpTbl.d5rt == ii),...
        'o', 'LineWidth', 1);
    e.Color = cols(cIdx(ii+1),:);
    e.MarkerFaceColor = cols(cIdx(ii+1),:);
    
    
    plot([1:nCohDisc]+xoff(ii), (grpTbl.mean_mean_accY(grpTbl.d5rt == ii)), '-', 'LineWidth', 2, 'Color', cols(cIdx(ii+1),:));
    
    allY = [allY; e.YData(:)];
    
end
set(gca, 'TickDir', 'out', 'LineWidth', 1)
xlim([.5, nCohDisc+.5])

yticks(ylim)
xticks([]);

yticks([ceil(min(allY)*100)/100, floor(max(allY)*100)/100])

title('Target \times RT')
xlabel('Target Coherence')

ylabel('relative accuracy')
yticks(ylim)
xticks([]);









% targ x RT (STATIC SIM)
grpTarg = grpstats(static_tbl, {'rt5', 'targ7'});
figM = nexttile; hold on;
for ii = 1:5
    plot((grpTarg.mean_acc(grpTarg.rt5 == ii)),...
        '-o', 'LineWidth', 2, 'Color', cols(cIdx(ii+1),:), 'MarkerFaceColor', cols(cIdx(ii+1),:));
end
set(gca, 'TickDir', 'out', 'LineWidth', 1);
title('static sim')
yticks(ylim)
xticks([]);
xlabel('target coherence')
ylabel('relative accuracy')
xlim([.5, nCohDisc+.5])




% allign axes
% lims = [min([figP.YAxis(1).Limits, figM.YAxis(1).Limits]), max([figP.YAxis(1).Limits, figM.YAxis(1).Limits])]
lims = [min([figP.YAxis(1).Limits, figM.YAxis(1).Limits]), .96]


figP.YAxis.Limits = lims; figP.YAxis.TickValues = lims;
figM.YAxis.Limits = lims; figM.YAxis.TickValues = lims; 











% ============== targ dynamics RT ==============
grpTbl = grpstats(rt_tbl, {'pt', 'acc', 'targ'}, {'mean', @skewness, @kurtosis});
grpTbl = grpstats(grpTbl, {'acc', 'targ'}, {'mean','sem'});




rt_cIdx = [180, 100];

xoff = [-.15, .15];

figP=nexttile([1,1]); hold on;
allY = [];
for ii = 0:1
    
    
        e=errorbar([1:nCohDisc]  + xoff(ii+1), (exp(grpTbl.mean_mean_rt(grpTbl.acc == ii)) + mean(all_ndt)), exp(grpTbl.mean_mean_rt(grpTbl.acc == ii)).*grpTbl.sem_mean_crt(grpTbl.acc == ii),...
            'o', 'LineWidth', 1);
    
%     e=errorbar([1:nCohDisc]  + xoff(ii+1), grpTbl.mean_skewness_rawRt(grpTbl.acc == ii), grpTbl.sem_skewness_rawRt(grpTbl.acc == ii),...
%         'o', 'LineWidth', 1);
    
    e.Color = cols(rt_cIdx(ii+1),:);
    e.MarkerFaceColor = cols(rt_cIdx(ii+1),:);
    
    
    plot([1:nCohDisc] + xoff(ii+1), (exp(grpTbl.mean_mean_rtY(grpTbl.acc == ii))+ mean(all_ndt)), '-', 'LineWidth', 2, 'Color', cols(rt_cIdx(ii+1),:));

    allY = [allY; e.YData(:)];
    
    
end
set(gca, 'TickDir', 'out', 'LineWidth', 1)
xlim([.5, nCohDisc+.5])
xticks([]);
title('Target \times Accuracy')

yticks([ceil(min(allY)*100)/100, floor(max(allY)*100)/100])

ylabel('RT')
yticks(ylim)
xticks([]);







% targ x Acc (STATIC SIM)
grpTarg = grpstats(static_tbl, {'acc', 'targ7'});
figM=nexttile; hold on;
for ii = 0:1
    plot((grpTarg.mean_rt(grpTarg.acc == ii)),...
        '-o', 'LineWidth', 2, 'Color', cols(rt_cIdx(ii+1),:), 'MarkerFaceColor', cols(rt_cIdx(ii+1),:));
end
set(gca, 'TickDir', 'out', 'LineWidth', 1);
title('static sim')
yticks(ylim)
xticks([]);
xlabel('target coherence')
ylabel('RT')
xlim([.5, nCohDisc+.5])



% allign axes
lims = [min([figP.YAxis(1).Limits, figM.YAxis(1).Limits]), max([figP.YAxis(1).Limits, figM.YAxis(1).Limits])]

figP.YAxis.Limits = lims; figP.YAxis.TickValues = lims;
figM.YAxis.Limits = lims; figM.YAxis.TickValues = lims; 





% ============ Distractor









% ============== dist dynamics acc ==============
grpTbl = grpstats(acc_tbl, {'pt', 'd5rt', 'dist'});
grpTbl = grpstats(grpTbl, {'d5rt', 'dist'}, {'mean','sem'});

figP=nexttile([1,1]); hold on;

xoff = linspace(-.2,.2,5);

allY = [];
for ii = 1:5
    
    e=errorbar([1:nCohDisc]+xoff(ii), center(grpTbl.mean_mean_acc(grpTbl.d5rt == ii)), grpTbl.sem_mean_cacc(grpTbl.d5rt == ii),...
        'o', 'LineWidth', 1);
    e.Color = cols(cIdx(ii+1),:);
    e.MarkerFaceColor = cols(cIdx(ii+1),:);
    
    plot([1:nCohDisc]+xoff(ii), center(grpTbl.mean_mean_accY(grpTbl.d5rt == ii)), '-', 'LineWidth', 2, 'Color', cols(cIdx(ii+1),:));
    
    allY = [allY; e.YData(:)];
    
end
set(gca, 'TickDir', 'out', 'LineWidth', 1)
xlim([.5, nCohDisc+.5])
xticks([]);

% legend('Fastest RT','','','','Slowest RT', 'Location', 'southeast')
ylabel('Relative Accuracy')
title('Distractor \times RT')
xlabel('distractor congruence')

yticks(ylim)
xticks([]);




% dist X RT  (STATIC SIM)
grpDist = grpstats(static_tbl, {'rt5', 'dist7'});
figM=nexttile; hold on;
for ii = 1:5
    plot(center(grpDist.mean_acc(grpDist.rt5 == ii)),...
        '-o', 'LineWidth', 2, 'Color', cols(cIdx(ii+1),:), 'MarkerFaceColor', cols(cIdx(ii+1),:));
end

set(gca, 'TickDir', 'out', 'LineWidth', 1);
title('static sim')
yticks(ylim)
xticks([]);
xlabel('distractor congruence')
ylabel('relative accuracy')
xlim([.5, nCohDisc+.5])



% allign axes
lims = [min([figP.YAxis(1).Limits, figM.YAxis(1).Limits]), max([figP.YAxis(1).Limits, figM.YAxis(1).Limits])]

figP.YAxis.Limits = lims; figP.YAxis.TickValues = lims;
figM.YAxis.Limits = lims; figM.YAxis.TickValues = lims; 








% ============== dist dynamics RT ==============
grpTbl = grpstats(rt_tbl, {'pt', 'acc', 'dist'});
grpTbl = grpstats(grpTbl, {'acc', 'dist'}, {'mean','sem'});

xoff = [-.15, .15];


figP=nexttile([1,1]); hold on;
allY = [];
for ii = 0:1
    
    
    e=errorbar([1:nCohDisc] + xoff(ii+1), (exp(grpTbl.mean_mean_rt(grpTbl.acc == ii)) + mean(all_ndt)), exp(grpTbl.mean_mean_rt(grpTbl.acc == ii)).*grpTbl.sem_mean_crt(grpTbl.acc == ii),...
        'o', 'LineWidth', 1);
    e.Color = cols(rt_cIdx(ii+1),:);
    e.MarkerFaceColor = cols(rt_cIdx(ii+1),:);
    
    
    plot([1:nCohDisc] + xoff(ii+1), (exp(grpTbl.mean_mean_rtY(grpTbl.acc == ii))+ mean(all_ndt)), '-', 'LineWidth', 2, 'Color', cols(rt_cIdx(ii+1),:));
    
    allY = [allY; e.YData(:)];
    
end
set(gca, 'TickDir', 'out', 'LineWidth', 1)
xlim([.5, nCohDisc+.5])
xticks([]);
xlabel('Distractor Congruence')

yticks([ceil(min(allY)*100)/100, floor(max(allY)*100)/100])

ylabel('RT')
title('Distractor \times Accuracy')
yticks(ylim)
xticks([]);






% dist x accuracy (STATIC SIM)
grpDist = grpstats(static_tbl, {'acc', 'dist7'});
figM=nexttile; hold on;
for ii = 0:1
    plot((grpDist.mean_rt(grpDist.acc == ii)),...
        '-o', 'LineWidth', 2, 'Color', cols(rt_cIdx(ii+1),:), 'MarkerFaceColor', cols(rt_cIdx(ii+1),:));
end

set(gca, 'TickDir', 'out', 'LineWidth', 1);
title('static sim')
yticks(ylim)
xticks([]);
xlabel('distractor congruence')
ylabel('RT')
xlim([.5, nCohDisc+.5])




% allign axes
lims = [min([figP.YAxis(1).Limits, figM.YAxis(1).Limits]), max([figP.YAxis(1).Limits, figM.YAxis(1).Limits])]

figP.YAxis.Limits = lims; figP.YAxis.TickValues = lims;
figM.YAxis.Limits = lims; figM.YAxis.TickValues = lims; 





%% plot dynamics trace ============================



load('batlow.mat');
cols = batlow;

% ===== simulation parameters  =====
threshType = 'race';    % which threshold {race, diff}
doRelu = 1;             % force positive accumulators?

nChoices = 2;           % how many accumulators?
nDt    = 2000;          % max number of timesteps
dt = .001;              % simulation dt

% ===== task parameters =====

nTarg = 11;
nDist = 11;
targRange = linspace(.6,.95,nTarg);
distRange = linspace(-1,1,nDist);

[targ, dist] = meshgrid(targRange, distRange);
targ = tanh(1.5*targ(:))./tanh(1.5);
ntarg = tanh(1.5*(1-targ(:)))./tanh(1.5);
dist = tanh(1.35*dist(:))./tanh(1.35);

nConds = length(targ);





% ========= run COLOR simulation




% decision dynamics

decay = 0;          % accumulator-specific decay
cross = 0;          % accumulator competition
ffinh = 1;          % feedforward inhibition
sd    = 1;          % diffusion noise
lapse = 0;          % lapse rate


A = cross*ones(nChoices) + (decay - cross)*eye(nChoices)
B = [...
    1, -ffinh;...
    -ffinh, 1;...
    ]


% gain dynamics
targ0       = 1;
targAim     = 6;
targGain    = 3;
 

dist0       = 1.25;
distAim     = 0;
distGain    = 4.5;

salSD       = 0.1;

 % non-decision time
t0 = .300;         

% decision threshold
thresh = max(.01, linspace(1.5, .1, nDt))';         



sddt = sd*sqrt(dt);
salSDdt = salSD*sqrt(dt);


% pre allocation




% run sims  ======================

clear meanRT Accuracy static_allRT static_allChoice

tic

x = zeros(nChoices,1);
all_x = nan(nDt,nChoices,1);
[all_targSal, all_distSal] = deal(nan(nDt,1));

%     targ0       = 2;
%     targAim     = 5;
%     targGain    = 1.5;
%
%     dist0       = 1;
%     distAim     = 0;
%     distGain    = 1.5;

targSal = targ0;
distSal = dist0;


cc = 121

I = [...
    targ(cc),  abs(dist(cc))*(dist(cc)>0);...
    ntarg(cc), abs(dist(cc))*(dist(cc)<0);...
    ];



for tt = 1:nDt
    
    % update gain
    targSal = targSal + targGain*(targAim - targSal)*dt + randn*salSDdt;
    distSal = distSal + distGain*(distAim - distSal)*dt + randn*salSDdt;
    
    % update model
    all_x(tt,:,:) = x;
    x = x + A*(x*dt) + B*(I*[targSal;distSal]*dt + randn(size(x))*sddt);
    
    % threshold
    if doRelu
        x = max(x,0);
    end
    
    
    
    all_targSal(tt) = targSal;
    all_distSal(tt) = distSal;
    all_x(tt,:,:) = x;
    
    % break if all accumulators reach threshold
    switch threshType
        
        case 'race'
            
            if ~any(max(x)<thresh(tt))
                break;
            end
            
        case 'diff'
            
            if ~any(abs(diff(x))<thresh(tt))
                break;
            end
            
    end
    
end



% keyboard;

finsel = isfinite(all_x(:,1));

dts = linspace(t0, sum(finsel)*dt + t0, sum(finsel));

figure('units','inch','position',[0,0,8,8]); hold on;
plot(dts, all_x(finsel,1), '-g', 'LineWidth', 2.5);
plot(dts, all_x(finsel,2), '-c', 'LineWidth', 2.5);
plot(dts, thresh(finsel), '--k', 'LineWidth', 1);
title('decision dynamics')
set(gca, 'TickDir', 'out', 'LineWidth', 1);





% === plot attractor dynamics

targSal = targ0;
distSal = dist0;



for tt = 1:1000
    
    targSal(tt+1) = targSal(tt) + targGain*(targAim - targSal(tt))*dt + randn*salSDdt;
    distSal(tt+1) = distSal(tt) + distGain*(distAim - distSal(tt))*dt + randn*salSDdt;
    
end





figure('units','inch','position',[0,0,16,8]); hold on;
tiledlayout(1,2, 'TileSpacing', 'compact', 'Padding', 'compact');


finsel = isfinite(targSal);
fin_finsel = find(finsel, 1, 'last');
dts = linspace(t0, sum(finsel)*dt + t0, sum(finsel));

nexttile; hold on;
plot(dts, targSal(finsel), '-g', 'LineWidth', 2.5);
plot(dts, distSal(finsel), '-c', 'LineWidth', 2.5);
yline(targAim, '--g', 'LineWidth', 1)
yline(distAim, '--c', 'LineWidth', 1)
title('salience dynamics')
set(gca, 'TickDir', 'out', 'LineWidth', 1);
ylabel('time')
ylabel('coherence gain')





    nexttile; hold on;

ndiscdt = 12;

discDt = round(linspace(1, length(targSal/3), ndiscdt));
coldtidx = round(linspace(1, length(cols), ndiscdt));

for ii = 1:ndiscdt
    
    plot(targSal(discDt(ii)), distSal(discDt(ii)), 'ok', 'LineWidth', 1, 'Color', cols(coldtidx(ii),:), 'MarkerFaceColor', cols(coldtidx(ii),:), 'MarkerSize', 10);
    
end




[xx,yy] = meshgrid(linspace(targ0-1,targAim+1,9), linspace(distAim-1, dist0+.75,9));
dTarg = targGain*(targAim - xx)*dt; 
dDist = distGain*(distAim - yy)*dt;
quiver(xx,yy,dTarg,dDist, 'k', 'LineWidth', .5);


xlabel('target gain')
ylabel('distractor gain')

title('attractor dynamics')
set(gca, 'TickDir', 'out', 'LineWidth', 1);
xlim([targ0-1, targAim+1])
ylim([distAim-1, dist0+.75])
xticks(xlim)
yticks(ylim)
yline(0, '-k', 'LineWidth', 1);









%%  FIT LCA within  


nReals = 10000;


% ========= run LCA simulation (run FFIc code block first)


targ0       = 3.5;
targAim     = 3.5;

dist0       = .4;
distAim     = .4;

t0 = .180;          % non-decision time
thresh = max(.01, linspace(2.5, -1, nDt))';         % decision threshold

lca_cross = -.60;
lca_decay = -.10;
lca_ffih  =    0;


a = lca_cross*ones(nChoices) + (lca_decay - lca_cross)*eye(nChoices);
sddt = sd*sqrt(dt);
salSDdt = salSD*sqrt(dt);


% run sims  ======================

clear stdRT meanRT Accuracy lca_allRT lca_allChoice

tic
parfor cc = 1:nConds
    
    x = zeros(nChoices,nReals);
    all_x = nan(nDt,nChoices,nReals);
    
    %     targ0       = 2;
    %     targAim     = 5;
    %     targGain    = 1.5;
    %
    %     dist0       = 1;
    %     distAim     = 0;
    %     distGain    = 1.5;
    
    targSal = targ0;
    distSal = dist0;
    
    for tt = 1:nDt
        
        targSal = targSal + targGain*(targAim - targSal)*dt + randn*salSDdt;
        distSal = distSal + distGain*(distAim - distSal)*dt + randn*salSDdt;
        
        
        I = [...
            targSal*targ(cc)  + distSal*abs(dist(cc))*(dist(cc)>0) - lca_ffih*( targSal*ntarg(cc) + distSal*abs(dist(cc))*(dist(cc)<0) );...
            targSal*ntarg(cc) + distSal*abs(dist(cc))*(dist(cc)<0) - lca_ffih*( targSal*targ(cc)  + distSal*abs(dist(cc))*(dist(cc)>0) ),...
            ];
        
        
        % update model
        all_x(tt,:,:) = x;
        x = x + (a*x + I)*dt + randn(size(x))*sddt;
        
        % threshold
        if doRelu
            x = max(x,0);
        end
        
        % break if all accumulators reach threshold
        switch threshType
            
            case 'race'
                
                if ~any(max(x)<thresh(tt))
                    break;
                end
                
            case 'diff'
                
                if ~any(abs(diff(x))<thresh(tt))
                    break;
                end
                
        end
        
    end
    
   [RT, Choice] = deal(nan(nReals, 1));
    
    for rr = 1:nReals
        try
            switch threshType
                case 'race'
                    RT(rr)      = find(max(all_x(:,:,rr),[],2)>=thresh, 1);
                case 'diff'
                    RT(rr)      = find(abs(diff(all_x(:,:,rr),[],2))>=thresh, 1);
            end
            
            [~,Choice(rr)]  = max(all_x(RT(rr),:,rr));
        catch
            RT(rr) = nan;
            Choice(rr) = 0;
        end
    end
    
    skewRT_cor(cc) = skewness(RT(Choice == 1 & isfinite(RT))*dt + t0);
    skewRT_err(cc) = skewness(RT(Choice == 2 & isfinite(RT))*dt + t0);
    stdRT(cc) = nanstd(RT(Choice == 1)*dt + t0);
    meanRT(cc) = nanmean(RT(Choice == 1)*dt + t0);
    Accuracy(cc) = nanmean(Choice == 1);
    
    lca_allRT{cc} = RT;
    lca_allChoice{cc} = Choice;
    
end
toc


% get mean
rStdRTLCA = reshape(stdRT,[nTarg, nDist]);
rRTLCA = reshape(meanRT,[nTarg, nDist]);
rAccLCA = reshape(Accuracy,[nTarg, nDist]);




% make table
[t_dist, t_targ, t_acc, t_rt] = deal([]);
for cc = 1:nConds
    
    t_acc  = [t_acc; lca_allChoice{cc}==1];
    t_rt   = [t_rt; lca_allRT{cc}*dt + t0];
    t_targ = [t_targ; ones(nReals,1)*targ(cc)];
    t_dist = [t_dist; ones(nReals,1)*dist(cc)];
    
end



lca_tbl = table;
lca_tbl.acc = t_acc;
lca_tbl.rt = t_rt;
lca_tbl.lrt = log(t_rt);
lca_tbl.crt = t_rt - nanmean(t_rt);
lca_tbl.rt5 = discretize(tiedrank(t_rt),5);
lca_tbl.targ = t_targ;
lca_tbl.targ7 = discretize(tiedrank(t_targ),7);
lca_tbl.dist = t_dist;
lca_tbl.dist7 = discretize(tiedrank(t_dist),7);



% regression
rt_mdl = fitlm(static_tbl, 'lrt ~ targ*dist', 'Exclude', static_tbl.acc == 0)
acc_mdl = fitglm(static_tbl, 'acc ~ targ*dist', 'Distribution', 'binomial')
dyn_mdl = fitglm(static_tbl, 'acc ~ (targ + dist)*crt', 'Distribution', 'binomial')
















load('batlow.mat');
cols = batlow;
cIdx = round(linspace(1, length(cols), 7));

rtacc_cols = {[36, 123, 160]/255, [242, 95, 92]/255};






mdl = wn_col;


[rt, crt, rt_pred, rt_targ, rt_dist, rt_res, rt_pt,...
    acc, cacc, acc_pred, acc_targ, acc_dist, acc_res, acc_pt,...
    all_ndt,...
    ] = deal([]);

pt_count = 1;

for mm = 1:length(col)
    
    
    data = mdl{mm}.data;
    results = mdl{mm}.results;
    
    param = results.E';
    ndt = ...
        results.link{ismember(results.paramName, 'rt_{ndt}')} ...
        (param(:, ismember(results.paramName, 'rt_{ndt}')));
    
    all_ndt = [all_ndt; ndt];
    
    for tt = 1:size(data,2)
        
        %                 if mm == 1 && tt == size(data,2)
        %                     continue
        %                 end
        
        % simulation parameters
        [~, out]=eval([results.likfun, '(param(tt,:), data(tt), [], [], 0, [])']);
        
        
        prt         = log(out.rt.rt(out.rt.rt_sel) - ndt(tt));
        prt(out.rt.rt(out.rt.rt_sel) <= ndt(tt)) = nan;
        rt          = [rt; prt];
        crt          = [crt; prt - nanmean(prt)];
        
        
        rt_pred     = [rt_pred; out.rt.rt_yhat];
        
        rt_targ     = [rt_targ; out.data.d11Targ0(out.rt.rt_sel)];
        rt_dist     = [rt_dist; 12-out.data.d11Dist0(out.rt.rt_sel)];
        %         rt_targ    = [rt_targ; out.data.d11Targ0(out.rt.rt_sel).*out.data.sgnResp(out.rt.rt_sel)];
        %         rt_dist    = [rt_dist; (12-out.data.d11Dist0(out.rt.rt_sel)).*out.data.sgnResp(out.rt.rt_sel)];
        
        rt_res      = [rt_res; prt - out.rt.rt_yhat];
        
        rt_pt       = [rt_pt; out.rt.rt_X(:,1).*pt_count];
        
        
        %         acc         = [acc; out.acc.acc(out.acc.acc_sel)];
        %         acc_pred    = [acc_pred; out.acc.acc_yhat];
        %         acc_targ    = [acc_targ; out.data.d11Targ0(out.acc.acc_sel).*out.data.sgnResp(out.acc.acc_sel)];
        %         acc_dist    = [acc_dist; (12-out.data.d11Dist0(out.acc.acc_sel)).*out.data.sgnResp(out.acc.acc_sel)];
        %         acc_pt      = [acc_pt; out.acc.acc_X(:,1).*pt_count];
        
        acc         = [acc; out.data.acc0(out.acc.acc_sel)];
        cacc        = [cacc; center(out.data.acc0(out.acc.acc_sel))];
        
        yhat        = [out.acc.acc_yhat];
        yhat(out.data.corrResp(out.acc.acc_sel)==0) = 1-yhat(out.data.corrResp(out.acc.acc_sel)==0);
        acc_pred    = [acc_pred; yhat];
        
        
        acc_targ    = [acc_targ; out.data.d11Targ0(out.acc.acc_sel)];
        acc_dist    = [acc_dist; (12-out.data.d11Dist0(out.acc.acc_sel))];
        acc_pt      = [acc_pt; out.acc.acc_X(:,1).*pt_count];
        
        pt_count = pt_count+1;
        
    end
    
    
end


cols = {[36, 123, 160]/255, [242, 95, 92]/255};



rt_tbl = array2table([rt, crt, rt_pred, rt_targ, rt_dist, rt_res, rt_pt ], 'VariableNames', ...
    {'rt', 'crt', 'rtY', 'targAll', 'distAll', 'res', 'pt'});

acc_tbl = array2table([acc, cacc, acc_pred, acc_targ, acc_dist, acc_pt], 'VariableNames', ...
    {'acc', 'cacc', 'accY', 'targAll', 'distAll', 'pt'});










% ===== MOTION


mdl = wn_mot;
% mdl = exp23;


[rt, crt, d5rt,rawRt, rt_pred, rt_targ, rt_dist, rt_acc, rt_pt,...
    acc, cacc, acc_pred, acc_targ, acc_dist, acc_rt, acc_pt,...
    all_ndtMot, rt_targ_all, rt_dist_all, acc_targ_all, acc_dist_all, ...
    ] = deal([]);

prevCor = [];

pt_count = 1;

nCohDisc = 7;


for mm = 1:3
    
    
    data = mdl{mm}.data;
    results = mdl{mm}.results;
    
    param = results.E';
    ndt = ...
        results.link{ismember(results.paramName, 'rt_{ndt}')} ...
        (param(:, ismember(results.paramName, 'rt_{ndt}')));
    
    all_ndtMot = [all_ndtMot; ndt];
    
    for tt = 1:size(data,2)
        
        %                 if mm == 1 && tt == size(data,2)
        %                     continue
        %                 end
        
        % simulation parameters
        [~, out]=eval([results.likfun, '(transpose(param(tt,:)), data(tt), [], [], 0, [])']);
        
        
        % rt
        prt         = log(out.rt.rt(out.rt.rt_sel) - ndt(tt));
        prt(out.rt.rt(out.rt.rt_sel) <= ndt(tt)) = nan;
        rt          = [rt; prt];
        crt          = [crt; prt - nanmean(prt)];
        rawRt       = [rawRt; out.rt.rt(out.rt.rt_sel)];
        
        
        rt_pred     = [rt_pred; out.rt.rt_yhat];
        
        rt_targ     = [acc_targ; discretize(tiedrank(out.data.targ0(out.rt.rt_sel)),nCohDisc)];
        rt_dist     = [acc_dist; discretize(tiedrank(out.data.dist0(out.rt.rt_sel)),nCohDisc)];
        
        
        rt_targ_all     = [rt_targ_all; out.data.d11Targ0(out.rt.rt_sel)];
        rt_dist_all     = [rt_dist_all; 12-out.data.d11Dist0(out.rt.rt_sel)];
        
        
        rt_acc      = [rt_acc; out.data.acc0(out.rt.rt_sel)];
        
        rt_pt       = [rt_pt; out.rt.rt_X(:,1).*pt_count];
        
        
        % acc
        d5rt          = [d5rt; discretize(tiedrank(out.data.rt0(out.acc.acc_sel)),5)];
        
        acc         = [acc; out.data.acc0(out.acc.acc_sel)];
        cacc        = [cacc; center(out.data.acc0(out.acc.acc_sel))];
        
        yhat        = [out.acc.acc_yhat];
        yhat(out.data.corrResp(out.acc.acc_sel)==0) = 1-yhat(out.data.corrResp(out.acc.acc_sel)==0);
        acc_pred    = [acc_pred; yhat];
        
        
        acc_targ     = [acc_targ; discretize(tiedrank(out.data.targ0(out.acc.acc_sel)),nCohDisc)];
        acc_dist     = [acc_dist; discretize(tiedrank(out.data.dist0(out.acc.acc_sel)),nCohDisc)];
        
        acc_targ_all     = [acc_targ_all; out.data.d11Targ0(out.acc.acc_sel)];
        acc_dist_all     = [acc_dist_all; 12-out.data.d11Dist0(out.acc.acc_sel)];
        
        acc_rt          = [acc_rt; out.data.rt0(out.acc.acc_sel)];
        acc_pt      = [acc_pt; out.acc.acc_X(:,1).*pt_count];
        
        pt_count = pt_count+1;
        
    end
    
    
end



rtMot_tbl = array2table([rt, crt, rt_pred, rt_targ_all, rt_dist_all, rt_acc, rt_pt ], 'VariableNames', ...
    {'rt', 'crt', 'rtY',  'targAll', 'distAll', 'acc', 'pt'});

accMot_tbl = array2table([acc, cacc, acc_pred, acc_rt, d5rt, acc_targ_all, acc_dist_all, acc_pt], 'VariableNames', ...
    {'acc', 'cacc', 'accY', 'rt', 'd5rt',  'targAll', 'distAll', 'pt'});










% acc_tbl = array2table([acc, acc_pred, rt, d5rt, acc_targ, acc_dist, acc_pt], 'VariableNames', ...
%     {'acc', 'accY','rt', 'd5rt', 'targ', 'dist' 'pt'});


npt = [size(col{1}.data,2), size(col{2}.data,2), size(col{3}.data,2)];
nptSel = [ones(npt(1),1)*1; ones(npt(2),1)*2; ones(npt(3),1)*3];








% ===== start plot
figure('units','inch','position',[0,0,2*7.5,7.5]); hold on;
tiledlayout(2,2, 'TileSpacing', 'compact', 'Padding', 'compact');






% =================== TARGET


nexttile([1,1]); hold on;

% ============== targ acc ==============
grpTbl = grpstats(acc_tbl, {'pt','targAll'});
grpTbl = grpstats(grpTbl, {'targAll'}, {'mean','sem'});

yyaxis('left');

% plot([1:length(grpTbl.mean_mean_acc)],(grpTbl.mean_mean_acc),'ob', 'LineWidth', .5);
ep = errorbar(1.15:11.15, grpTbl.mean_mean_acc, grpTbl.sem_mean_cacc,...
    'o', 'color', rtacc_cols{1}, 'MarkerFaceColor', rtacc_cols{1}, 'LineWidth', 1);

plot(1.15:11.15, grpTbl.mean_mean_accY, '-', 'color', rtacc_cols{1}, 'LineWidth', 2);
set(gca, 'TickDir', 'out', 'LineWidth', 1)
yticks([floor(min(ep.YData)*100)/100, ceil(max(ep.YData)*100)/100])

colTargAccY = ylim;
yticks(colTargAccY)
xticks([]);

% ============== targ rt ==============
grpTbl = grpstats(rt_tbl, {'pt','targAll'}, {'mean'});
grpTbl = grpstats(grpTbl, {'targAll'}, {'mean','sem'});

yyaxis('right');

% plot(exp(grpTbl.mean_mean_rt) + mean(all_ndt), 'ob', 'LineWidth', .5);
ep=errorbar(.85:10.85, exp(grpTbl.mean_mean_rt) + mean(all_ndt), exp(grpTbl.mean_mean_rt).*(grpTbl.sem_mean_crt),...
    'o', 'color', rtacc_cols{2}, 'MarkerFaceColor', rtacc_cols{2}, 'LineWidth', 1);

plot(.85:10.85, exp(grpTbl.mean_mean_rtY) + mean(all_ndt), '-', 'color', rtacc_cols{2}, 'LineWidth', 2);

set(gca, 'TickDir', 'out', 'LineWidth', 1)
yticks([floor(min(ep.YData)*100)/100, ceil(max(ep.YData)*100)/100])


title('Target (Attend-Color)')

xlim([.5 11.5])
xticks([])
xticklabels({'low coh', 'high coh'})

xlabel('target coherence')


colTargRTY = ylim;
yticks(colTargRTY)
xticks([]);





% plot target LCA color =========

% 
% % target
% nexttile([1,1]); hold on;
% 
% yyaxis('left');
% plot(.85:10.85, mean(rAccLCA,1), '-o', 'color', rtacc_cols{1}, 'MarkerFaceColor', rtacc_cols{1}, 'LineWidth', 2);
% ylabel('accuracy')
% xticks([]);
% yticks([]);
% xlim([.5 11.5])
% 
% yyaxis('right');
% plot(.85:10.85, mean(rRTLCA,1), '-o', 'color', rtacc_cols{2}, 'MarkerFaceColor', rtacc_cols{2}, 'LineWidth', 2);
% ylabel('RT')
% 
% set(gca, 'LineWidth', 1, 'TickDir', 'out');
% title('LCA sim')
% xlabel('target coherence')
% xticks([]);
% yticks([]);
% xlim([.5 11.5])









% plot target LCA color =========



% target
nexttile([1,1]); hold on;

yyaxis('left');
plot(.85:10.85, mean(rAccLCA,1), '-o', 'color', rtacc_cols{1}, 'MarkerFaceColor', rtacc_cols{1}, 'LineWidth', 2);
ylabel('accuracy')

% ylim(colTargAccY);
yticks(ylim)
xticks([]);



yyaxis('right');
plot(.85:10.85, mean(rRTLCA,1), '-o', 'color', rtacc_cols{2}, 'MarkerFaceColor', rtacc_cols{2}, 'LineWidth', 2);
ylabel('RT')

set(gca, 'LineWidth', 1, 'TickDir', 'out');
title('DDM sim')
xlabel('target coherence')

% ylim(colTargRTY);
yticks(ylim)
xticks([]);













% distractor


nexttile([1,1]); hold on;

% ============== dist acc ==============
yyaxis('left');

grpTbl = grpstats(acc_tbl, {'pt','distAll'});
grpTbl = grpstats(grpTbl, {'distAll'}, {'mean','sem'});

% plot([1:length(grpTbl.mean_mean_acc)], (grpTbl.mean_mean_acc),'ok', 'LineWidth', .5);
ep=errorbar(1.15:11.15, grpTbl.mean_mean_acc, grpTbl.sem_mean_cacc,...
    'o', 'color', rtacc_cols{1}, 'MarkerFaceColor', rtacc_cols{1}, 'LineWidth', 1);
plot(1.15:11.15, grpTbl.mean_mean_accY, '-', 'color', rtacc_cols{1}, 'LineWidth', 2);
set(gca, 'TickDir', 'out', 'LineWidth', 1)
yticks([floor(min(ep.YData)*100)/100, ceil(max(ep.YData)*100)/100])
xticks([])


motDistAccY = ylim;
yticks(motDistAccY)
xticks([]);


% ============== dist rt ==============
yyaxis('right');

grpTbl = grpstats(rt_tbl, {'pt','distAll'}, {'mean'});
grpTbl = grpstats(grpTbl, {'distAll'}, {'mean','sem'});



% plot(exp(grpTbl.mean_mean_rt) + mean(all_ndt), 'ok', 'LineWidth', .5);
ep = errorbar(.85:10.85, exp(grpTbl.mean_mean_rt) + mean(all_ndt), exp(grpTbl.mean_mean_rt).*(grpTbl.sem_mean_crt),...
    'o', 'color', rtacc_cols{2}, 'MarkerFaceColor', rtacc_cols{2}, 'LineWidth', 1);
plot(.85:10.85, exp(grpTbl.mean_mean_rtY)+ mean(all_ndt), '-', 'color', rtacc_cols{2}, 'LineWidth', 2);
set(gca, 'TickDir', 'out', 'LineWidth', 1)
title('Distractor (Attend-Color)')
% ylim(targRange);

yticks([floor(min(ep.YData)*100)/100, ceil(max(ep.YData)*100)/100])
xlim([.5 11.5])
xticks([])
xlabel('distractor congruence')


motDistRTY = ylim;
yticks(ylim)
xticks([]);













% plot distractor LCA color =========



nexttile([1,1]); hold on;

yyaxis('left');
plot(.85:10.85, mean(rAccLCA,2), '-o', 'color', rtacc_cols{1}, 'MarkerFaceColor', rtacc_cols{1}, 'LineWidth', 2);
ylabel('accuracy')


% ylim(motDistAccY);
yticks(ylim)
xticks([]);


xlim([.5 11.5])

yyaxis('right');
plot(.85:10.85, mean(rRTLCA,2), '-o', 'color', rtacc_cols{2}, 'MarkerFaceColor', rtacc_cols{2}, 'LineWidth', 2);
ylabel('RT')

set(gca, 'LineWidth', 1, 'TickDir', 'out');
title('DDM sim')
xlabel('distractor congruence')


% ylim(motDistRTY);
yticks(ylim)
xticks([]);









%% FFIc Rew & Adapt LOAD

clear adaptDyn rewDyn 

adaptDyn{1} = load('2021-04-27_17-25_exp1_adaptDualDyn.mat');
adaptDyn{2} = load('2021-02-02_06-22_exp2_adaptDualDyn.mat');
adaptDyn{3} = load('2021-02-03_02-03_exp3_adaptDualDyn.mat');

rewDyn{3} = load('2021-04-06_20-56_exp3_rewDualDyn.mat');


disp('loaded')




%%  ~~~~ Updated FFIc dynamic [Rew & Adapt] ~~~~~~~

nReals = 10000;          % how many simulated trials?



% ===== simulation parameters  =====
threshType = 'race';    % which threshold {race, diff}
doRelu = 1;             % force positive accumulators?

nChoices = 2;           % how many accumulators?
nDt    = 2000;          % max number of timesteps
dt = .001;              % simulation dt

% ===== task parameters =====

nTarg = 11;
nDist = 11;
targRange = linspace(.6,.95,nTarg);
distRange = linspace(-1,1,nDist);

[targ, dist] = meshgrid(targRange, distRange);
targ = tanh(1.5*targ(:))./tanh(1.5);
ntarg = tanh(1.5*(1-targ(:)))./tanh(1.5);
dist = tanh(1.35*dist(:))./tanh(1.35);

nConds = length(targ);








% decision dynamics

decay = 0;          % accumulator-specific decay
cross = 0;          % accumulator competition
ffinh = 1;          % feedforward inhibition
sd    = 1;          % diffusion noise
lapse = 0;          % lapse rate


A = cross*ones(nChoices) + (decay - cross)*eye(nChoices)
B = [...
    1, -ffinh;...
    -ffinh, 1;...
    ]


% gain dynamics
targAim     = 6;
targGain    = 3;
 

distAim     = 0;
distGain    = 4.5;

salSD       = 0.1;

 % non-decision time
t0 = .350;         

% decision threshold
thresh = max(.01, linspace(1.5, .1, nDt))';         




sddt = sd*sqrt(dt);
salSDdt = salSD*sqrt(dt);







% ========= run ADAPT LOW simulation

targ0       = 1;
dist0       = 1.5;




% run sims  ======================

clear stdRT meanRT Accuracy allRT allChoice rSkewRT_cor rSkewRT_err

tic
parfor cc = 1:nConds
    
    x = zeros(nChoices,nReals);
    all_x = nan(nDt,nChoices,nReals);
    
    %     targ0       = 2;
    %     targAim     = 5;
    %     targGain    = 1.5;
    %
    %     dist0       = 1;
    %     distAim     = 0;
    %     distGain    = 1.5;
    
    targSal = targ0;
    distSal = dist0;
    

    I = [...
        targ(cc),  abs(dist(cc))*(dist(cc)>0);...
        ntarg(cc), abs(dist(cc))*(dist(cc)<0);...
        ];
    
    
    for tt = 1:nDt
        
        % update gain
        targSal = targSal + targGain*(targAim - targSal)*dt + randn*salSDdt;
        distSal = distSal + distGain*(distAim - distSal)*dt + randn*salSDdt;
        
        % update model
        all_x(tt,:,:) = x;
        x = x + A*(x*dt) + B*(I*[targSal;distSal]*dt + randn(size(x))*sddt);
        
        % threshold
        if doRelu
            x = max(x,0);
        end
        
        % break if all accumulators reach threshold
        switch threshType
            
            case 'race'
                
                if ~any(max(x)<thresh(tt))
                    break;
                end
                
            case 'diff'
                
                if ~any(abs(diff(x))<thresh(tt))
                    break;
                end
                
        end
        
    end
    
    [RT, Choice] = deal(nan(nReals, 1));
    
    for rr = 1:nReals
        try
            switch threshType
                case 'race'
                    RT(rr)      = find(max(all_x(:,:,rr),[],2) >= thresh, 1);
                case 'diff'
                    RT(rr)      = find(abs(diff(all_x(:,:,rr),[],2)) >= thresh, 1);
            end
            
            [~,Choice(rr)]  = max(all_x(RT(rr),:,rr));
            
        catch
            RT(rr) = nan;
            Choice(rr) = 0;
        end
    end
    
    skewRT_cor(cc) = skewness(RT(Choice == 1 & isfinite(RT))*dt + t0);
    skewRT_err(cc) = skewness(RT(Choice == 2 & isfinite(RT))*dt + t0);
    stdRT(cc) = nanstd(RT(Choice == 1)*dt + t0);
    meanRT(cc) = nanmean(RT(Choice == 1)*dt + t0);
    Accuracy(cc) = nanmean(Choice == 1);
    
    allRT{cc} = RT;
    allChoice{cc} = Choice;
    
end
toc


% get mean
rSkewRT_cor = reshape(skewRT_cor,[nTarg, nDist]);
rSkewRT_err = reshape(skewRT_err,[nTarg, nDist]);
rStdRT = reshape(stdRT,[nTarg, nDist]);
rRT = reshape(meanRT,[nTarg, nDist]);
rAcc = reshape(Accuracy,[nTarg, nDist]);



% make table
[t_dist, t_targ, t_acc, t_rt] = deal([]);
for cc = 1:nConds
    
    t_acc  = [t_acc; allChoice{cc}==1];
    t_rt   = [t_rt; allRT{cc}*dt + t0];
    t_targ = [t_targ; ones(nReals,1)*targ(cc)];
    t_dist = [t_dist; ones(nReals,1)*dist(cc)];
    
end

mean(isfinite(t_rt))


tbl = table;
tbl.acc = t_acc;
tbl.rt = t_rt;


% ntr = height(tbl);
% nLapse = round(ntr*.05);
% lapseTr = randsample(1:ntr, nLapse);
% 
% tbl.acc(lapseTr) = round(rand(nLapse,1));
% tbl.rt(lapseTr) = unifrnd(t0, 2, nLapse,1);


tbl.lrt = log(t_rt);
tbl.crt = t_rt - nanmean(t_rt);
tbl.rt5 = discretize(tiedrank(t_rt),5);
tbl.targ = t_targ;
tbl.targ5 = discretize(tiedrank(t_targ),5);
tbl.dist = t_dist;
tbl.dist5 = discretize(tiedrank(t_dist),5);











% ADAPT HIGH  ==============================================================
targ0       = 1;
dist0       = 1.5 - 1.25;


% run sims  ======================

clear stdRT meanRT Accuracy allRT allChoice rSkewRT_cor rSkewRT_err

tic
parfor cc = 1:nConds
    
    x = zeros(nChoices,nReals);
    all_x = nan(nDt,nChoices,nReals);
    
    %     targ0       = 2;
    %     targAim     = 5;
    %     targGain    = 1.5;
    %
    %     dist0       = 1;
    %     distAim     = 0;
    %     distGain    = 1.5;
    
    targSal = targ0;
    distSal = dist0;
    

    I = [...
        targ(cc),  abs(dist(cc))*(dist(cc)>0);...
        ntarg(cc), abs(dist(cc))*(dist(cc)<0);...
        ];
    
    
    for tt = 1:nDt
        
        % update gain
        targSal = targSal + targGain*(targAim - targSal)*dt + randn*salSDdt;
        distSal = distSal + distGain*(distAim - distSal)*dt + randn*salSDdt;
        
        % update model
        all_x(tt,:,:) = x;
        x = x + A*(x*dt) + B*(I*[targSal;distSal]*dt + randn(size(x))*sddt);
        
        % threshold
        if doRelu
            x = max(x,0);
        end
        
        % break if all accumulators reach threshold
        switch threshType
            
            case 'race'
                
                if ~any(max(x)<thresh(tt))
                    break;
                end
                
            case 'diff'
                
                if ~any(abs(diff(x))<thresh(tt))
                    break;
                end
                
        end
        
    end
    
    [RT, Choice] = deal(nan(nReals, 1));
    
    for rr = 1:nReals
        try
            switch threshType
                case 'race'
                    RT(rr)      = find(max(all_x(:,:,rr),[],2) >= thresh, 1);
                case 'diff'
                    RT(rr)      = find(abs(diff(all_x(:,:,rr),[],2)) >= thresh, 1);
            end
            
            [~,Choice(rr)]  = max(all_x(RT(rr),:,rr));
            
        catch
            RT(rr) = nan;
            Choice(rr) = 0;
        end
    end
    
    skewRT_cor(cc) = skewness(RT(Choice == 1 & isfinite(RT))*dt + t0);
    skewRT_err(cc) = skewness(RT(Choice == 2 & isfinite(RT))*dt + t0);
    stdRT(cc) = nanstd(RT(Choice == 1)*dt + t0);
    meanRT(cc) = nanmean(RT(Choice == 1)*dt + t0);
    Accuracy(cc) = nanmean(Choice == 1);
    
    allRT{cc} = RT;
    allChoice{cc} = Choice;
    
end
toc


% get mean
rSkewRT_cor = reshape(skewRT_cor,[nTarg, nDist]);
rSkewRT_err = reshape(skewRT_err,[nTarg, nDist]);
rStdRT = reshape(stdRT,[nTarg, nDist]);
rRT = reshape(meanRT,[nTarg, nDist]);
rAcc = reshape(Accuracy,[nTarg, nDist]);



% make table
[t_dist, t_targ, t_acc, t_rt] = deal([]);
for cc = 1:nConds
    
    t_acc  = [t_acc; allChoice{cc}==1];
    t_rt   = [t_rt; allRT{cc}*dt + t0];
    t_targ = [t_targ; ones(nReals,1)*targ(cc)];
    t_dist = [t_dist; ones(nReals,1)*dist(cc)];
    
end

mean(isfinite(t_rt))


adaptSim_tbl = table;
adaptSim_tbl.acc = t_acc;
adaptSim_tbl.rt = t_rt;


% ntr = height(adaptSim_tbl);
% nLapse = round(ntr*.05);
% lapseTr = randsample(1:ntr, nLapse);
% 
% adaptSim_tbl.acc(lapseTr) = round(rand(nLapse,1));
% adaptSim_tbl.rt(lapseTr) = unifrnd(t0, 2, nLapse,1);


adaptSim_tbl.lrt = log(t_rt);
adaptSim_tbl.crt = t_rt - nanmean(t_rt);
adaptSim_tbl.rt5 = discretize(tiedrank(t_rt),5);
adaptSim_tbl.targ = t_targ;
adaptSim_tbl.targ5 = discretize(tiedrank(t_targ),5);
adaptSim_tbl.dist = t_dist;
adaptSim_tbl.dist5 = discretize(tiedrank(t_dist),5);










% REW LOW  ==============================================================
targ0       = .90;
dist0       = 1.25;

thresh = max(.01, linspace(1, -.2, nDt))';         % decision threshold



clear stdRT meanRT Accuracy allRT allChoice rSkewRT_cor rSkewRT_err

tic
parfor cc = 1:nConds
    
    x = zeros(nChoices,nReals);
    all_x = nan(nDt,nChoices,nReals);
    
    %     targ0       = 2;
    %     targAim     = 5;
    %     targGain    = 1.5;
    %
    %     dist0       = 1;
    %     distAim     = 0;
    %     distGain    = 1.5;
    
    targSal = targ0;
    distSal = dist0;
    

    I = [...
        targ(cc),  abs(dist(cc))*(dist(cc)>0);...
        ntarg(cc), abs(dist(cc))*(dist(cc)<0);...
        ];
    
    
    for tt = 1:nDt
        
        % update gain
        targSal = targSal + targGain*(targAim - targSal)*dt + randn*salSDdt;
        distSal = distSal + distGain*(distAim - distSal)*dt + randn*salSDdt;
        
        % update model
        all_x(tt,:,:) = x;
        x = x + A*(x*dt) + B*(I*[targSal;distSal]*dt + randn(size(x))*sddt);
        
        % threshold
        if doRelu
            x = max(x,0);
        end
        
        % break if all accumulators reach threshold
        switch threshType
            
            case 'race'
                
                if ~any(max(x)<thresh(tt))
                    break;
                end
                
            case 'diff'
                
                if ~any(abs(diff(x))<thresh(tt))
                    break;
                end
                
        end
        
    end
    
    [RT, Choice] = deal(nan(nReals, 1));
    
    for rr = 1:nReals
        try
            switch threshType
                case 'race'
                    RT(rr)      = find(max(all_x(:,:,rr),[],2) >= thresh, 1);
                case 'diff'
                    RT(rr)      = find(abs(diff(all_x(:,:,rr),[],2)) >= thresh, 1);
            end
            
            [~,Choice(rr)]  = max(all_x(RT(rr),:,rr));
            
        catch
            RT(rr) = nan;
            Choice(rr) = 0;
        end
    end
    
    skewRT_cor(cc) = skewness(RT(Choice == 1 & isfinite(RT))*dt + t0);
    skewRT_err(cc) = skewness(RT(Choice == 2 & isfinite(RT))*dt + t0);
    stdRT(cc) = nanstd(RT(Choice == 1)*dt + t0);
    meanRT(cc) = nanmean(RT(Choice == 1)*dt + t0);
    Accuracy(cc) = nanmean(Choice == 1);
    
    allRT{cc} = RT;
    allChoice{cc} = Choice;
    
end
toc


% get mean
rSkewRT_cor = reshape(skewRT_cor,[nTarg, nDist]);
rSkewRT_err = reshape(skewRT_err,[nTarg, nDist]);
rStdRT = reshape(stdRT,[nTarg, nDist]);
rRT = reshape(meanRT,[nTarg, nDist]);
rAcc = reshape(Accuracy,[nTarg, nDist]);



% make table
[t_dist, t_targ, t_acc, t_rt] = deal([]);
for cc = 1:nConds
    
    t_acc  = [t_acc; allChoice{cc}==1];
    t_rt   = [t_rt; allRT{cc}*dt + t0];
    t_targ = [t_targ; ones(nReals,1)*targ(cc)];
    t_dist = [t_dist; ones(nReals,1)*dist(cc)];
    
end

mean(isfinite(t_rt))


lowRew = table;
lowRew.acc = t_acc;
lowRew.rt = t_rt;


% ntr = height(lowRew);
% nLapse = round(ntr*.05);
% lapseTr = randsample(1:ntr, nLapse);
% 
% lowRew.acc(lapseTr) = round(rand(nLapse,1));
% lowRew.rt(lapseTr) = unifrnd(t0, 2, nLapse,1);


lowRew.lrt = log(t_rt);
lowRew.crt = t_rt - nanmean(t_rt);
lowRew.rt5 = discretize(tiedrank(t_rt),5);
lowRew.targ = t_targ;
lowRew.targ5 = discretize(tiedrank(t_targ),5);
lowRew.dist = t_dist;
lowRew.dist5 = discretize(tiedrank(t_dist),5);












% REW HIGH  ==============================================================
targ0       = 1.45;
dist0       = 1.25;

% targ0       = 1.5;
% targAim     = 1.5;



clear stdRT meanRT Accuracy allRT allChoice rSkewRT_cor rSkewRT_err

tic
parfor cc = 1:nConds
    
    x = zeros(nChoices,nReals);
    all_x = nan(nDt,nChoices,nReals);
    
    %     targ0       = 2;
    %     targAim     = 5;
    %     targGain    = 1.5;
    %
    %     dist0       = 1;
    %     distAim     = 0;
    %     distGain    = 1.5;
    
    targSal = targ0;
    distSal = dist0;
    

    I = [...
        targ(cc),  abs(dist(cc))*(dist(cc)>0);...
        ntarg(cc), abs(dist(cc))*(dist(cc)<0);...
        ];
    
    
    for tt = 1:nDt
        
        % update gain
        targSal = targSal + targGain*(targAim - targSal)*dt + randn*salSDdt;
        distSal = distSal + distGain*(distAim - distSal)*dt + randn*salSDdt;
        
        % update model
        all_x(tt,:,:) = x;
        x = x + A*(x*dt) + B*(I*[targSal;distSal]*dt + randn(size(x))*sddt);
        
        % threshold
        if doRelu
            x = max(x,0);
        end
        
        % break if all accumulators reach threshold
        switch threshType
            
            case 'race'
                
                if ~any(max(x)<thresh(tt))
                    break;
                end
                
            case 'diff'
                
                if ~any(abs(diff(x))<thresh(tt))
                    break;
                end
                
        end
        
    end
    
    [RT, Choice] = deal(nan(nReals, 1));
    
    for rr = 1:nReals
        try
            switch threshType
                case 'race'
                    RT(rr)      = find(max(all_x(:,:,rr),[],2) >= thresh, 1);
                case 'diff'
                    RT(rr)      = find(abs(diff(all_x(:,:,rr),[],2)) >= thresh, 1);
            end
            
            [~,Choice(rr)]  = max(all_x(RT(rr),:,rr));
            
        catch
            RT(rr) = nan;
            Choice(rr) = 0;
        end
    end
    
    skewRT_cor(cc) = skewness(RT(Choice == 1 & isfinite(RT))*dt + t0);
    skewRT_err(cc) = skewness(RT(Choice == 2 & isfinite(RT))*dt + t0);
    stdRT(cc) = nanstd(RT(Choice == 1)*dt + t0);
    meanRT(cc) = nanmean(RT(Choice == 1)*dt + t0);
    Accuracy(cc) = nanmean(Choice == 1);
    
    allRT{cc} = RT;
    allChoice{cc} = Choice;
    
end
toc


% get mean
rSkewRT_cor = reshape(skewRT_cor,[nTarg, nDist]);
rSkewRT_err = reshape(skewRT_err,[nTarg, nDist]);
rStdRT = reshape(stdRT,[nTarg, nDist]);
rRT = reshape(meanRT,[nTarg, nDist]);
rAcc = reshape(Accuracy,[nTarg, nDist]);



% make table
[t_dist, t_targ, t_acc, t_rt] = deal([]);
for cc = 1:nConds
    
    t_acc  = [t_acc; allChoice{cc}==1];
    t_rt   = [t_rt; allRT{cc}*dt + t0];
    t_targ = [t_targ; ones(nReals,1)*targ(cc)];
    t_dist = [t_dist; ones(nReals,1)*dist(cc)];
    
end

mean(isfinite(t_rt))


rewSim_tbl = table;
rewSim_tbl.acc = t_acc;
rewSim_tbl.rt = t_rt;




rewSim_tbl.lrt = log(t_rt);
rewSim_tbl.crt = t_rt - nanmean(t_rt);
rewSim_tbl.rt5 = discretize(tiedrank(t_rt),5);
rewSim_tbl.targ = t_targ;
rewSim_tbl.targ5 = discretize(tiedrank(t_targ),5);
rewSim_tbl.dist = t_dist;
rewSim_tbl.dist5 = discretize(tiedrank(t_dist),5);

















% ====== LOAD ADAPT ========



load('batlow.mat');
bat_cols = batlow;
bat_cIdx = round(linspace(1, length(bat_cols), 7));


mdl = adaptDyn;


[rt, crt, rt_acc, rt_pred, rt_targ, rt_dist, rt_pt,...
    rt_rew, acc_rew,...
    acc, cacc, acc_pred, acc_targ, acc_dist, acc_pt,...
    all_ndt, exp2_ndt,...
    rt_blk, acc_blk, acc_sBlk, acc_rt, d5rt acc_dist1,acc_crt, ...
    ] = deal([]);

pt_count = 1;

whichExps = 1:3;

nBlkDisc = 2;
nCohDisc = 3;



rtacc_cols = {[36, 123, 160]/255, [242, 95, 92]/255};


for mm = whichExps
    
    data = mdl{mm}.data;
    results = mdl{mm}.results;
    
    param = results.E';
    ndt = ...
        results.link{ismember(results.paramName, 'rt_{ndt}')} ...
        (param(:, ismember(results.paramName, 'rt_{ndt}')));
    
    if mm == 3
        all_ndt = [all_ndt; ndt];
    else
        exp2_ndt = [exp2_ndt; ndt];
    end
    
    for tt = 1:size(data,2)
        
        % simulation parameters
        [~, out]=eval([results.likfun, '(param(tt,:), data(tt), [], [], 0, [])']);
        
        
        prt         = log((out.rt.rt(out.rt.rt_sel) - ndt(tt)));
        prt(out.rt.rt(out.rt.rt_sel) <= ndt(tt)) = nan;
        rt          = [rt; prt];
        crt         = [crt; prt - nanmean(prt)];
        
        rt_acc      = [rt_acc; out.data.acc0(out.rt.rt_sel)];
        
        rt_pred     = [rt_pred; out.rt.rt_yhat];
        
        %                 rt_targ     = [rt_targ; discretize(tiedrank(out.data.targ0(out.rt.rt_sel)),nCohDisc)];
        %                 rt_dist     = [rt_dist; discretize(tiedrank(out.data.dist0(out.rt.rt_sel)),nCohDisc)];
        
        rt_targ     = [rt_targ; out.data.d11Targ0(out.rt.rt_sel)];
        rt_dist     = [rt_dist; -out.data.d11Dist0(out.rt.rt_sel)];
        
        if mm == 2
            rt_rew      = [rt_rew; out.data.d11Targ0(out.rt.rt_sel)*0];
        else
            rt_rew      = [rt_rew; grp2idx(out.data.rew(out.rt.rt_sel))];
        end
        
        rt_pt       = [rt_pt; out.rt.rt_X(:,1).*pt_count];
        
        
        
        %         acc         = [acc; out.acc.acc(out.acc.acc_sel)];
        %         acc_pred    = [acc_pred; out.acc.acc_yhat];
        %         acc_targ    = [acc_targ; out.data.d11Targ0(out.acc.acc_sel).*out.data.sgnResp(out.acc.acc_sel)];
        %         acc_dist    = [acc_dist; (12-out.data.d11Dist0(out.acc.acc_sel)).*out.data.sgnResp(out.acc.acc_sel)];
        %         acc_pt      = [acc_pt; out.acc.acc_X(:,1).*pt_count];
        
        acc         = [acc; out.data.acc0(out.acc.acc_sel)];
        cacc        = [cacc; center(out.data.acc0(out.acc.acc_sel))];
        
        yhat        = [out.acc.acc_yhat];
        yhat(out.data.corrResp(out.acc.acc_sel)==0) = 1-yhat(out.data.corrResp(out.acc.acc_sel)==0);
        acc_pred    = [acc_pred; yhat];
        
        acc_rt          = [acc_rt; out.data.rt0(out.acc.acc_sel)];
        acc_crt          = [acc_crt; out.acc.crt];

        
        
        condDRT    = out.data.rt0(out.acc.acc_sel);
        condDRT(grp2idx(out.data.rew(out.acc.acc_sel)) == 1)    = discretize(tiedrank(condDRT(grp2idx(out.data.rew(out.acc.acc_sel)) == 1)),5);
        condDRT(grp2idx(out.data.rew(out.acc.acc_sel)) == 2)    = discretize(tiedrank(condDRT(grp2idx(out.data.rew(out.acc.acc_sel)) == 2)),5);
        
        d5rt          = [d5rt; discretize(tiedrank(out.data.rt0(out.acc.acc_sel)),5)];
%         d5rt          = [d5rt; condDRT];

        
        %                 acc_targ    = [acc_targ; discretize(tiedrank(out.data.targ0(out.acc.acc_sel)), nCohDisc)];
        %                 acc_dist    = [acc_dist; discretize(tiedrank((out.data.dist0(out.acc.acc_sel))), nCohDisc)];
        
        acc_targ    = [acc_targ; discretize(tiedrank((out.data.targ0(out.acc.acc_sel))), 5)];
        acc_dist    = [acc_dist; discretize(tiedrank((out.data.dist0(out.acc.acc_sel))), 5)];
        
        acc_dist1    = [acc_dist1; discretize(tiedrank((out.data.dist1(out.acc.acc_sel))), nCohDisc)];

        
        if mm == 2
            acc_rew      = [acc_rew; out.data.d11Targ0(out.acc.acc_sel)*0];
        else
            acc_rew      = [acc_rew; grp2idx(out.data.rew(out.acc.acc_sel))];
        end
        
        acc_pt      = [acc_pt; out.acc.acc_X(:,1).*pt_count];
        pt_count = pt_count+1;
        
        
        
        
        
        
    end
    
    
end


dyn_rt_tbl = array2table([rt, rt_acc, crt, rt_pred, rt_targ, rt_dist, rt_rew, rt_pt ], 'VariableNames', ...
    {'rt', 'acc', 'crt', 'rtY', 'targ', 'dist', 'rew', 'pt'});

adapt_tbl = array2table([acc, cacc, acc_pred, acc_targ, acc_dist, acc_rt, acc_crt, d5rt, acc_rew, acc_dist1, acc_pt], 'VariableNames', ...
    {'acc', 'cacc', 'accY', 'targ', 'dist', 'rt', 'crt', 'drt', 'rew', 'dist1', 'pt'});











% ====== LOAD REW ========




load('batlow.mat');
bat_cols = batlow;
bat_cIdx = round(linspace(1, length(bat_cols), 7));


mdl = rewDyn;


[rt, crt, rt_acc, rt_pred, rt_targ, rt_dist, rt_pt,...
    rt_rew, acc_rew,...
    acc, cacc, acc_pred, acc_targ, acc_dist, acc_pt,...
    all_ndt, exp2_ndt,...
    rt_blk, acc_blk, acc_sBlk, acc_rt, acc_crt, d5rt ...
    ] = deal([]);

pt_count = 1;

whichExps = 3;

nBlkDisc = 2;
nCohDisc = 5;



rtacc_cols = {[36, 123, 160]/255, [242, 95, 92]/255};


for mm = whichExps
    
    data = mdl{mm}.data;
    results = mdl{mm}.results;
    
    param = results.E';
    ndt = ...
        results.link{ismember(results.paramName, 'rt_{ndt}')} ...
        (param(:, ismember(results.paramName, 'rt_{ndt}')));
    
    if mm == 3
        all_ndt = [all_ndt; ndt];
    else
        exp2_ndt = [exp2_ndt; ndt];
    end
    
    for tt = 1:size(data,2)
        
        % simulation parameters
        [~, out]=eval([results.likfun, '(param(tt,:), data(tt), [], [], 0, [])']);
        
        
        prt         = log((out.rt.rt(out.rt.rt_sel) - ndt(tt)));
        prt(out.rt.rt(out.rt.rt_sel) <= ndt(tt)) = nan;
        rt          = [rt; prt];
        crt         = [crt; prt - nanmean(prt)];
        
        rt_acc      = [rt_acc; out.data.acc0(out.rt.rt_sel)];
        
        rt_pred     = [rt_pred; out.rt.rt_yhat];
        
        %                 rt_targ     = [rt_targ; discretize(tiedrank(out.data.targ0(out.rt.rt_sel)),nCohDisc)];
        %                 rt_dist     = [rt_dist; discretize(tiedrank(out.data.dist0(out.rt.rt_sel)),nCohDisc)];
        
        rt_targ     = [rt_targ; out.data.d11Targ0(out.rt.rt_sel)];
        rt_dist     = [rt_dist; -out.data.d11Dist0(out.rt.rt_sel)];
        
        if mm == 2
            rt_rew      = [rt_rew; out.data.d11Targ0(out.rt.rt_sel)*0];
        else
            rt_rew      = [rt_rew; grp2idx(out.data.rew(out.rt.rt_sel))];
        end
        
        rt_pt       = [rt_pt; out.rt.rt_X(:,1).*pt_count];
        
        
        
        %         acc         = [acc; out.acc.acc(out.acc.acc_sel)];
        %         acc_pred    = [acc_pred; out.acc.acc_yhat];
        %         acc_targ    = [acc_targ; out.data.d11Targ0(out.acc.acc_sel).*out.data.sgnResp(out.acc.acc_sel)];
        %         acc_dist    = [acc_dist; (12-out.data.d11Dist0(out.acc.acc_sel)).*out.data.sgnResp(out.acc.acc_sel)];
        %         acc_pt      = [acc_pt; out.acc.acc_X(:,1).*pt_count];
        
        acc         = [acc; out.data.acc0(out.acc.acc_sel)];
        cacc        = [cacc; center(out.data.acc0(out.acc.acc_sel))];
        
        yhat        = [out.acc.acc_yhat];
        yhat(out.data.corrResp(out.acc.acc_sel)==0) = 1-yhat(out.data.corrResp(out.acc.acc_sel)==0);
        acc_pred    = [acc_pred; yhat];
        
        acc_rt          = [acc_rt; out.data.rt0(out.acc.acc_sel)];
        acc_crt         = [acc_crt; out.acc.crt];

        
        
        condDRT    = out.data.rt0(out.acc.acc_sel);
        condDRT(grp2idx(out.data.rew(out.acc.acc_sel)) == 1)    = discretize(tiedrank(condDRT(grp2idx(out.data.rew(out.acc.acc_sel)) == 1)),5);
        condDRT(grp2idx(out.data.rew(out.acc.acc_sel)) == 2)    = discretize(tiedrank(condDRT(grp2idx(out.data.rew(out.acc.acc_sel)) == 2)),5);
        
        d5rt          = [d5rt; discretize(tiedrank(out.data.rt0(out.acc.acc_sel)),5)];
%         d5rt          = [d5rt; condDRT];

        
        %                 acc_targ    = [acc_targ; discretize(tiedrank(out.data.targ0(out.acc.acc_sel)), nCohDisc)];
        %                 acc_dist    = [acc_dist; discretize(tiedrank((out.data.dist0(out.acc.acc_sel))), nCohDisc)];
        
        acc_targ    = [acc_targ; out.data.d5Targ0(out.acc.acc_sel)];
        acc_dist    = [acc_dist; -out.data.d5Dist0(out.acc.acc_sel)];
        
        if mm == 2
            acc_rew      = [acc_rew; out.data.d11Targ0(out.acc.acc_sel)*0];
        else
            acc_rew      = [acc_rew; grp2idx(out.data.rew(out.acc.acc_sel))];
        end
        
        acc_pt      = [acc_pt; out.acc.acc_X(:,1).*pt_count];
        pt_count = pt_count+1;
        
        
        
        
        
        
    end
    
    
end


dyn_rt_tbl = array2table([rt, rt_acc, crt, rt_pred, rt_targ, rt_dist, rt_rew, rt_pt ], 'VariableNames', ...
    {'rt', 'acc', 'crt', 'rtY', 'targ', 'dist', 'rew', 'pt'});

rew_tbl = array2table([acc, cacc, acc_pred, acc_targ, acc_dist, acc_rt, acc_crt, d5rt, acc_rew, acc_pt], 'VariableNames', ...
    {'acc', 'cacc', 'accY', 'targ', 'dist', 'rt', 'crt', 'drt', 'rew', 'pt'});



















% ~~~~~ plot BOTH

% plot
figure('units','inch','position',[0,0,2*7.5,2*2]); hold on;
tiledlayout(1,4, 'TileSpacing', 'none', 'Padding', 'none');











% ============== plot Adapt Interaction Behav ==============



grpTbl = grpstats(adapt_tbl, {'pt', 'dist', 'drt' 'dist1'});
grpTbl = grpstats(grpTbl, {'dist', 'drt', 'dist1'}, {'mean','sem'});


nexttile([1,1]); hold on;
% yline(0, '--k', 'LineWidth', 1)

xxoff = linspace(-.1,.1,5);
for ii = [1:5]
    
%     plot([1:5]+xxoff(ii), (grpTbl.mean_mean_acc(grpTbl.drt == ii & grpTbl.dist1==3))-(grpTbl.mean_mean_acc(grpTbl.drt == ii & grpTbl.dist1==1)),...
%         'o', 'LineWidth', 1, 'color', bat_cols(bat_cIdx(ii+1),:),  'MarkerFaceColor', bat_cols(bat_cIdx(ii+1),:));
%     
%     plot([1:5]+xxoff(ii), (grpTbl.mean_mean_accY(grpTbl.drt == ii & grpTbl.dist1==3))-(grpTbl.mean_mean_accY(grpTbl.drt == ii & grpTbl.dist1==1)),...
%         '-', 'LineWidth', 2, 'color', bat_cols(bat_cIdx(ii+1),:));
    
    
     plot([1:5]+xxoff(ii), (grpTbl.mean_mean_acc(grpTbl.drt == ii & grpTbl.dist1==3)),...
        '--', 'LineWidth', 2, 'color', bat_cols(bat_cIdx(ii+1),:));
    
     plot([1:5]+xxoff(ii), (grpTbl.mean_mean_acc(grpTbl.drt == ii & grpTbl.dist1==1)),...
        '-', 'LineWidth', 2, 'color', bat_cols(bat_cIdx(ii+1),:));
    
    
end

% ylim([-.2, .1])
yticks(ylim)
set(gca, 'TickDir', 'out', 'LineWidth', 1);
xlim([.5, 5.5])
title('adaptation (BEHAV)')
ylabel('prev Cong - Inc (accuracy)')
xlabel('distractor congruence')
% legend(pltp, {'fast','','med','','slow'}, 'Location', 'northeast')







% ============== plot Adapt Interaction SIM ==============

adaptAll_tbl = [tbl; adaptSim_tbl];
adaptAll_tbl.rt5 = discretize(tiedrank(adaptAll_tbl.rt), 5);
adaptAll_tbl.adapt = [zeros(height(tbl),1); ones(height(adaptSim_tbl),1)];

grpTbl = grpstats(adaptAll_tbl, {'dist5', 'rt5', 'adapt'});


nexttile([1,1]); hold on;
% yline(0, '--k', 'LineWidth', 1)

xxoff = linspace(-.1,.1,5);
for ii = [1:5]
    
%     plot([1:5]+xxoff(ii), grpTbl.mean_acc(grpTbl.rt5 == ii & grpTbl.adapt == 0) - grpTbl.mean_acc(grpTbl.rt5 == ii & grpTbl.adapt == 1),...
%         '-o', 'LineWidth', 1.5, 'color', bat_cols(bat_cIdx(ii+1),:),  'MarkerFaceColor', bat_cols(bat_cIdx(ii+1),:));
    
    
        plot([1:5]+xxoff(ii), grpTbl.mean_acc(grpTbl.rt5 == ii & grpTbl.adapt == 1),...
        '-o', 'LineWidth', 1.5, 'color', bat_cols(bat_cIdx(ii+1),:),  'MarkerFaceColor', bat_cols(bat_cIdx(ii+1),:));
    
    
        plot([1:5]+xxoff(ii), grpTbl.mean_acc(grpTbl.rt5 == ii & grpTbl.adapt == 0),...
        '--o', 'LineWidth', 1.5, 'color', bat_cols(bat_cIdx(ii+1),:),  'MarkerFaceColor', bat_cols(bat_cIdx(ii+1),:));
   
     
end


% ylim([-.2, .1])
yticks(ylim)
set(gca, 'TickDir', 'out', 'LineWidth', 1);
xlim([.5, 5.5])
title('adaptation (SIM)')
ylabel('prev Cong - Inc (accuracy)')
xlabel('distractor congruence')















% ============== plot REW interaction ==============




grpTbl = grpstats(rew_tbl, {'pt', 'targ', 'drt' 'rew'});
grpTbl = grpstats(grpTbl, {'targ', 'drt', 'rew'}, {'mean','sem'});


nexttile([1,1]); hold on;
% yline(0, '--k', 'LineWidth', 1)

xxoff = linspace(-.1,.1,5);
for ii = [1:5]
    
%     plot([1:5]+xxoff(ii), grpTbl.mean_mean_acc(grpTbl.drt == ii & grpTbl.rew == 2)-grpTbl.mean_mean_acc(grpTbl.drt == ii & grpTbl.rew == 1),...
%         'o', 'LineWidth', 1, 'color', bat_cols(bat_cIdx(ii+1),:),  'MarkerFaceColor', bat_cols(bat_cIdx(ii+1),:));
%     
%     pltp(ii) = plot([1:5]+xxoff(ii), grpTbl.mean_mean_accY(grpTbl.drt == ii & grpTbl.rew == 2)-grpTbl.mean_mean_accY(grpTbl.drt == ii & grpTbl.rew == 1),...
%         '-', 'LineWidth', 2, 'color', bat_cols(bat_cIdx(ii+1),:));


    pltp(ii) = plot([1:5]+xxoff(ii), grpTbl.mean_mean_acc(grpTbl.drt == ii & grpTbl.rew == 2),...
        '-', 'LineWidth', 2, 'color', bat_cols(bat_cIdx(ii+1),:));

    pltp(ii) = plot([1:5]+xxoff(ii), grpTbl.mean_mean_acc(grpTbl.drt == ii & grpTbl.rew == 1),...
        '--', 'LineWidth', 2, 'color', bat_cols(bat_cIdx(ii+1),:));
 
          
end

% ylim([-.02, .09])

yticks(ylim)
set(gca, 'TickDir', 'out', 'LineWidth', 1);
xlim([.5, 5.5])
title('reward (BEHAV)')
ylabel('Hi - No REW (accuracy)')
xlabel('target coherence')









% ============== plot REW Interaction SIM ==============


rewAll_tbl = [lowRew; rewSim_tbl];
rewAll_tbl.rt5 = discretize(tiedrank(rewAll_tbl.rt), 5);
rewAll_tbl.rew = [zeros(height(lowRew),1); ones(height(adaptSim_tbl),1)];

grpTbl = grpstats(rewAll_tbl, {'dist5', 'rt5', 'rew'});


nexttile([1,1]); hold on;
% yline(0, '--k', 'LineWidth', 1)

xxoff = linspace(-.1,.1,5);
for ii = [1:5]
    
%     plot([1:5]+xxoff(ii), grpTbl.mean_acc(grpTbl.rt5 == ii & grpTbl.rew == 1) - grpTbl.mean_acc(grpTbl.rt5 == ii & grpTbl.rew == 0),...
%         '-o', 'LineWidth', 1.5, 'color', bat_cols(bat_cIdx(ii+1),:),  'MarkerFaceColor', bat_cols(bat_cIdx(ii+1),:));
    
    
    
        plot([1:5]+xxoff(ii), grpTbl.mean_acc(grpTbl.rt5 == ii & grpTbl.rew == 1),...
        '-o', 'LineWidth', 1.5, 'color', bat_cols(bat_cIdx(ii+1),:),  'MarkerFaceColor', bat_cols(bat_cIdx(ii+1),:));
    
    
        plot([1:5]+xxoff(ii), grpTbl.mean_acc(grpTbl.rt5 == ii & grpTbl.rew == 0),...
        '--o', 'LineWidth', 1.5, 'color', bat_cols(bat_cIdx(ii+1),:),  'MarkerFaceColor', bat_cols(bat_cIdx(ii+1),:));
    
    
end

% ylim([-.02, .09])
yticks(ylim)
set(gca, 'TickDir', 'out', 'LineWidth', 1);
xlim([.5, 5.5])
title('reward (SIM)')
ylabel('Hi - No REW (accuracy)')
xlabel('target coherence')












% ~~~~~ plot CONTRAST

% plot
figure('units','inch','position',[0,0,2*7.5,2*2]); hold on;
tiledlayout(1,4, 'TileSpacing', 'none', 'Padding', 'none');











% ============== plot Adapt Interaction Behav ==============



grpTbl = grpstats(adapt_tbl, {'pt', 'dist', 'drt' 'dist1'});
grpTbl = grpstats(grpTbl, {'dist', 'drt', 'dist1'}, {'mean','sem'});


nexttile([1,1]); hold on;
yline(0, '--k', 'LineWidth', 1)

xxoff = linspace(-.1,.1,5);
for ii = [1:5]
    
    plot([1:5]+xxoff(ii), (grpTbl.mean_mean_acc(grpTbl.drt == ii & grpTbl.dist1==3))-(grpTbl.mean_mean_acc(grpTbl.drt == ii & grpTbl.dist1==1)),...
        'o', 'LineWidth', 1, 'color', bat_cols(bat_cIdx(ii+1),:),  'MarkerFaceColor', bat_cols(bat_cIdx(ii+1),:));
    
    plot([1:5]+xxoff(ii), (grpTbl.mean_mean_accY(grpTbl.drt == ii & grpTbl.dist1==3))-(grpTbl.mean_mean_accY(grpTbl.drt == ii & grpTbl.dist1==1)),...
        '-', 'LineWidth', 2, 'color', bat_cols(bat_cIdx(ii+1),:));
    

end

% ylim([-.2, .1])
yticks(ylim)
set(gca, 'TickDir', 'out', 'LineWidth', 1);
xlim([.5, 5.5])
title('adaptation (BEHAV)')
ylabel('prev Cong - Inc (accuracy)')
xlabel('distractor congruence')
% legend(pltp, {'fast','','med','','slow'}, 'Location', 'northeast')







% ============== plot Adapt Interaction SIM ==============

adaptAll_tbl = [tbl; adaptSim_tbl];
adaptAll_tbl.rt5 = discretize(tiedrank(adaptAll_tbl.rt), 5);
adaptAll_tbl.adapt = [zeros(height(tbl),1); ones(height(adaptSim_tbl),1)];

grpTbl = grpstats(adaptAll_tbl, {'dist5', 'rt5', 'adapt'});


nexttile([1,1]); hold on;
yline(0, '--k', 'LineWidth', 1)

xxoff = linspace(-.1,.1,5);
for ii = [1:5]
    
    plot([1:5]+xxoff(ii), grpTbl.mean_acc(grpTbl.rt5 == ii & grpTbl.adapt == 0) - grpTbl.mean_acc(grpTbl.rt5 == ii & grpTbl.adapt == 1),...
        '-o', 'LineWidth', 1.5, 'color', bat_cols(bat_cIdx(ii+1),:),  'MarkerFaceColor', bat_cols(bat_cIdx(ii+1),:));
    
    
     
end


ylim([-.2, .1])
yticks(ylim)
set(gca, 'TickDir', 'out', 'LineWidth', 1);
xlim([.5, 5.5])
title('adaptation (SIM)')
ylabel('prev Cong - Inc (accuracy)')
xlabel('distractor congruence')















% ============== plot REW interaction ==============




grpTbl = grpstats(rew_tbl, {'pt', 'targ', 'drt' 'rew'});
grpTbl = grpstats(grpTbl, {'targ', 'drt', 'rew'}, {'mean','sem'});


nexttile([1,1]); hold on;
yline(0, '--k', 'LineWidth', 1)

xxoff = linspace(-.1,.1,5);
for ii = [1:5]
    
    plot([1:5]+xxoff(ii), grpTbl.mean_mean_acc(grpTbl.drt == ii & grpTbl.rew == 2)-grpTbl.mean_mean_acc(grpTbl.drt == ii & grpTbl.rew == 1),...
        'o', 'LineWidth', 1, 'color', bat_cols(bat_cIdx(ii+1),:),  'MarkerFaceColor', bat_cols(bat_cIdx(ii+1),:));
    
    pltp(ii) = plot([1:5]+xxoff(ii), grpTbl.mean_mean_accY(grpTbl.drt == ii & grpTbl.rew == 2)-grpTbl.mean_mean_accY(grpTbl.drt == ii & grpTbl.rew == 1),...
        '-', 'LineWidth', 2, 'color', bat_cols(bat_cIdx(ii+1),:));


end

ylim([-.02, .09])

yticks(ylim)
set(gca, 'TickDir', 'out', 'LineWidth', 1);
xlim([.5, 5.5])
title('reward (BEHAV)')
ylabel('Hi - No REW (accuracy)')
xlabel('target coherence')









% ============== plot REW Interaction SIM ==============


rewAll_tbl = [lowRew; rewSim_tbl];
rewAll_tbl.rt5 = discretize(tiedrank(rewAll_tbl.rt), 5);
rewAll_tbl.rew = [zeros(height(tbl),1); ones(height(adaptSim_tbl),1)];

grpTbl = grpstats(rewAll_tbl, {'dist5', 'rt5', 'rew'});


nexttile([1,1]); hold on;
yline(0, '--k', 'LineWidth', 1)

xxoff = linspace(-.1,.1,5);
for ii = [1:5]
    
    plot([1:5]+xxoff(ii), grpTbl.mean_acc(grpTbl.rt5 == ii & grpTbl.rew == 1) - grpTbl.mean_acc(grpTbl.rt5 == ii & grpTbl.rew == 0),...
        '-o', 'LineWidth', 1.5, 'color', bat_cols(bat_cIdx(ii+1),:),  'MarkerFaceColor', bat_cols(bat_cIdx(ii+1),:));
    
    
    
end

ylim([-.02, .09])
yticks(ylim)
set(gca, 'TickDir', 'out', 'LineWidth', 1);
xlim([.5, 5.5])
title('reward (SIM)')
ylabel('Hi - No REW (accuracy)')
xlabel('target coherence')


















% === plot attractor dynamics =============================================

figure; hold on;
nexttile; hold on;





% == orig
targ0       = 1.5;
dist0       = 1.25;

targSal = targ0;
distSal = dist0;



for tt = 1:800
    
    targSal(tt+1) = targSal(tt) + targGain*(targAim - targSal(tt))*dt + randn*salSDdt;
    distSal(tt+1) = distSal(tt) + distGain*(distAim - distSal(tt))*dt + randn*salSDdt;
    
end




finsel = isfinite(targSal);
fin_finsel = find(finsel, 1, 'last');
dts = linspace(t0, sum(finsel)*dt + t0, sum(finsel));

plot(dts, targSal(finsel), '-g', 'LineWidth', 2.5);
plot(dts, distSal(finsel), '-c', 'LineWidth', 2.5);
yline(targAim, '--g', 'LineWidth', 1)
yline(distAim, '--c', 'LineWidth', 1)
title('salience dynamics')
set(gca, 'TickDir', 'out', 'LineWidth', 1);
ylabel('time')
ylabel('coherence gain')







% == mod
targ0       = 1.5 + .5;
dist0       = 1.25 - 1;

targSal = targ0;
distSal = dist0;



for tt = 1:800
    
    targSal(tt+1) = targSal(tt) + targGain*(targAim - targSal(tt))*dt + randn*salSDdt;
    distSal(tt+1) = distSal(tt) + distGain*(distAim - distSal(tt))*dt + randn*salSDdt;
    
end




finsel = isfinite(targSal);
fin_finsel = find(finsel, 1, 'last');
dts = linspace(t0, sum(finsel)*dt + t0, sum(finsel));

plot(dts, targSal(finsel), '-k', 'LineWidth', 2.5);
plot(dts, distSal(finsel), '-k', 'LineWidth', 2.5);
yline(targAim, '--g', 'LineWidth', 1)
yline(distAim, '--c', 'LineWidth', 1)
title('salience dynamics')
set(gca, 'TickDir', 'out', 'LineWidth', 1);
ylabel('time')
ylabel('coherence gain')





