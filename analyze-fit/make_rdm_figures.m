%% === make RDM figures  ===
%
% run script from `analyze-fit` folder
%
%
% Harrison Ritz - Oct 2021




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






%% ====== TARGET & DISTRACTOR =============================================
clear col mot dyn


% tanh AR softplus
col{1} = load('2021-01-17_21-31_exp1_withinStatic.mat');
col{2} = load('2021-01-17_18-43_exp2_withinStatic.mat');
col{3} = load('2021-01-20_00-43_exp3_withinStatic.mat');

noX{1} = load('2021-01-17_21-31_exp1_withinStatic.mat');
noX{2} = load('2021-02-19_15-56_exp2_withinStaticNoX.mat');
noX{3} = load('2021-02-19_13-59_exp3_withinStaticNoX.mat');

% tanh AR softplus
mot{1} = load('2021-01-17_21-40_exp1_motionStatic.mat');
mot{2} = load('2021-01-17_22-28_exp2_motionStatic.mat');
mot{3} = load('2021-01-19_10-12_exp3_motionStatic.mat');

% tanh AR softplus
dyn{1} = load('2021-04-26_16-36_exp1_withinDualDyn.mat');
dyn{2} = load('2021-01-17_20-55_exp2_withinDualDyn.mat');
dyn{3} = load('2021-01-17_22-52_exp3_withinDualDyn.mat');






disp('');
disp('loaded');
disp('');



% save stats
clc;
disp_mdls = {'col', 'mot', 'dyn'};
mdlNames = {'attend-color', 'attend-motion', 'dynamics'};

for cc = 1:3
    
    for mm = 1:length(disp_mdls)
        
        % transform parameters
        x = eval([disp_mdls{mm},'{cc}.results.stats.alpha']);
        for xx = 1:length(x)
           x(xx) =  eval([disp_mdls{mm},'{cc}.data(1).link{xx}(x(xx))']);
        end
        
        
        % correct parameter names
        n = eval([disp_mdls{mm},'{cc}.results.paramName'])';
        for nn = 1:length(n)
            n(nn) = strrep(n(nn), 'acc', 'choice');
            n(nn) = strrep(n(nn), 'pwr', 'preGain');
            n(nn) = strrep(n(nn), 'targDist', 'targXdist');
            n(nn) = strrep(n(nn), 'crt', 'RT');
            n(nn) = strrep(n(nn), 'Crt', 'RT');
             n(nn) = strrep(n(nn), 'rt_{choice}', 'rt_{acc}');
            n(nn) = strrep(n(nn), 'rt_{choiceTarg}', 'rt_{accTarg}');
            n(nn) = strrep(n(nn), 'rt_{choiceDist}', 'rt_{accDist}');

        end
        %         disp(n);
        
        % == make table
        tbl = table;
        
        tbl.paramName           = n;
        tbl.df                  = eval(['repmat(length(',disp_mdls{mm},'{cc}.data) - length(', disp_mdls{mm},'{cc}.results.stats.tval), [ length(', disp_mdls{mm},'{cc}.results.stats.tval), 1])']);
        tbl.paramVal            = x;
        tbl.tVal                = eval([disp_mdls{mm},'{cc}.results.stats.tval']);
        tbl.pVal                = eval([disp_mdls{mm},'{cc}.results.stats.p']);
        tbl.cohenD              = eval([disp_mdls{mm},'{cc}.results.stats.alpha./sqrt(diag(', disp_mdls{mm}, '{cc}.results.stats.groupvar))']);
        
        
        
        
        
        
        fprintf('\n =============== exp %d %s =============== \n\n', cc, disp_mdls{mm});
        %         disp(tbl)
        
        writetable(tbl, './tables/withinTrial.xlsx', 'Sheet', sprintf('exp%d %s', cc, mdlNames{mm}))
        
    end
    
    %      fprintf('\n\n\n\n\n');
    
end



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

    combP1 = 2*normcdf(nansum(norminv(1-((realmin+pvals)./2)).*df, 2) ./ sqrt(isfinite(pvals)*(df.^2)'), 'upper');
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





for cc = 1:3
    
    [~, paramR] = cov2corr(col{cc}.results.stats.groupvar);
    
    ns = length(col{cc}.data);
    
    
    fprintf('\nexp %d (df = %f)', cc, ns-2);
    
    
    targ_r = paramR(ismember(col{cc}.results.paramName, 'rt_{targ}'), ismember(col{cc}.results.paramName, 'acc_{targ}'))
    targ_p = tcdf(-abs(targ_r * sqrt((ns-2)/(1-targ_r.^2))), ns-2)
    
    dist_r = paramR(ismember(col{cc}.results.paramName, 'rt_{dist}'), ismember(col{cc}.results.paramName, 'acc_{dist}'))
    dist_p = tcdf(-abs(dist_r * sqrt((ns-2)/(1-dist_r.^2))), ns-2)
    
    
end




% % ====== compare color & motion
% 
% for cc = 1:3
%     
%     npt = length(col{cc}.data);
%     
%     % acc targ
%     tcrt = ismember(col{cc}.results.paramName, 'acc_{targ}');
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
%     tcrt = ismember(col{cc}.results.paramName, 'acc_{dist}');
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
%     % rt targ
%     tcrt = ismember(col{cc}.results.paramName, 'rt_{targ}');
%     cv2 = col{cc}.results.stats.groupmeancovariance(tcrt,tcrt);
%     cv3 = mot{cc}.results.stats.groupmeancovariance(tcrt,tcrt);
%     
%     rtTarg_tval = (col{cc}.results.alpha(tcrt) - mot{cc}.results.alpha(tcrt))./...
%         sqrt(cv2 + cv3);
%     
%     rtTarg_df = (cv2 + cv3).^2 ./...
%         ((cv2.^2)./(npt-1) + (cv3.^2)./(npt-1));
%     
%     rtTarg_pval = 2*tcdf(abs(rtTarg_tval), rtTarg_df, 'upper');
%     
%     
%     % rt targ
%     tcrt = ismember(col{cc}.results.paramName, 'rt_{dist}');
%     cv2 = col{cc}.results.stats.groupmeancovariance(tcrt,tcrt);
%     cv3 = mot{cc}.results.stats.groupmeancovariance(tcrt,tcrt);
%     
%     rtDist_tval = (col{cc}.results.alpha(tcrt) - mot{cc}.results.alpha(tcrt))./...
%         sqrt(cv2 + cv3);
%     
%     rtDist_df = (cv2 + cv3).^2 ./...
%         ((cv2.^2)./(npt-1) + (cv3.^2)./(npt-1));
%     
%     rtDist_pval = 2*tcdf(abs(rtDist_tval), rtDist_df, 'upper');
%     
%     
%     fprintf('\n exp %d', cc');
%     fprintf('\n == mot-col targ acc: t(%.5g)=%.5g, p=%.5g',  accTarg_df, accTarg_tval, accTarg_pval)
%     fprintf('\n == mot-col dist acc: t(%.5g)=%.5g, p=%.5g',  accDist_df, accDist_tval, accDist_pval)
%     fprintf('\n == mot-col targ rt: t(%.5g)=%.5g, p=%.5g',  rtTarg_df, rtTarg_tval, rtTarg_pval)
%     fprintf('\n == mot-col dist rt: t(%.5g)=%.5g, p=%.5g\n\n',  rtDist_df, rtDist_tval, rtDist_pval)
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
    tcrt = ismember(col{cc}.results.paramName, 'acc_{dist}');
    mtcrt = ismember(mot{cc}.results.paramName, 'acc_{dist}');
    
    cv2 = col{cc}.results.stats.groupmeancovariance(tcrt,tcrt);
    cv3 = mot{cc}.results.stats.groupmeancovariance(mtcrt,mtcrt);
    
    accTarg_tval = (col{cc}.results.alpha(tcrt) - mot{cc}.results.alpha(mtcrt))./...
        sqrt(cv2 + cv3);
    
    accTarg_df = (col_idf*cv2 + mot_idf*cv3).^2 ./...
        (col_idf.*(col_idf*cv2).^2 + mot_idf.*(mot_idf*cv3).^2);
    
    accTarg_pval = 2*tcdf(abs(accTarg_tval), accTarg_df, 'upper');
    
    
    % acc targ
    tcrt = ismember(col{cc}.results.paramName, 'acc_{targ}');
    mtcrt = ismember(mot{cc}.results.paramName, 'acc_{targ}');
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
    tcrt = ismember(col{cc}.results.paramName, 'rt_{dist}');
    mtcrt = ismember(mot{cc}.results.paramName, 'rt_{dist}');
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
    tcrt = ismember(col{cc}.results.paramName, 'rt_{targ}');
    mtcrt = ismember(mot{cc}.results.paramName, 'rt_{targ}');
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
    fprintf('\n == col-mot dist acc: t(%.5g)=%.5g, d=%.5g, p=%.5g',  accTarg_df, accTarg_tval, accTarg_tval./sqrt(accTarg_df), accTarg_pval)
    fprintf('\n == col-mot targ acc: t(%.5g)=%.5g, d=%.5g, p=%.5g',  accDist_df, accDist_tval, accDist_tval./sqrt(accDist_df), accDist_pval)
    fprintf('\n == col-mot dist rt: t(%.5g)=%.5g, d=%.5g, p=%.5g',  rtTarg_df, rtTarg_tval, rtTarg_tval./sqrt(rtTarg_df), rtTarg_pval)
    fprintf('\n == col-mot targ rt: t(%.5g)=%.5g, d=%.5g, p=%.5g\n\n',  rtDist_df, rtDist_tval, rtDist_tval./sqrt(rtDist_df), rtDist_pval)
    
    if cc==1
        [accDist_tval, accDist_pval, accDist_df, rtDist_tval, rtDist_pval, rtDist_df] =deal(nan);
    end
    
    
    pvals(:,cc) = [accTarg_pval;  accDist_pval; rtTarg_pval; rtDist_pval];
    df(:,cc) = [accTarg_df; accDist_df; rtTarg_df; rtDist_df];
    
    sqrn(:,cc) = [cc_sqn; cc_sqn; cc_sqn; cc_sqn];
    sgnD(:,cc) = sign([accTarg_tval; accDist_tval; rtTarg_tval; rtDist_tval]);
    
    
end

sgnD = sgnD == sign(sgnD(:,end));


% combP = 2*normcdf(nansum(norminv(1-((realmin+pvals)./2)).*df, 2) ./ sqrt(nansum((df.^2),2)), 'upper')
combP = normcdf(nansum(sgnD.*-norminv(pvals).*sqn, 2) ./ sqrt(isfinite(pvals)*(sqn.^2)'), 'upper')









% display correlations

for ii = 2:3
    
    fprintf('\nexp %d\n', ii)
    
    [~,cv] = cov2corr(col{ii}.results.stats.groupvar);

    betweenFeatureCorr  =  [cv(ismember(col{ii}.data(1).paramName, 'rt_{targ}'), ismember(col{ii}.data(1).paramName, 'rt_{dist}')),...
        cv(ismember(col{ii}.data(1).paramName, 'acc_{targ}'), ismember(col{ii}.data(1).paramName, 'acc_{dist}'))]
    
    btFeatureP = tcdf(abs((betweenFeatureCorr*sqrt(df(ii))) ./ (sqrt(1-betweenFeatureCorr.^2))), df(ii), 'upper')*2
    
    
end


%% plot within-trial static


mdl = col;


[rt, crt, rt_pred, rt_targ, rt_dist, rt_targ5, rt_dist5, rt_res, rt_pt, rt_tr,...
    acc, cacc, acc_pred, acc_targ, acc_dist, acc_targ5, acc_dist5, acc_res, acc_pt, acc_tr,...
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
        
        rt_targ5     = [rt_targ5; out.data.d5Targ0(out.rt.rt_sel)];
        rt_dist5     = [rt_dist5; 6-out.data.d5Dist0(out.rt.rt_sel)];
        
        %         rt_targ    = [rt_targ; out.data.d11Targ0(out.rt.rt_sel).*out.data.sgnResp(out.rt.rt_sel)];
        %         rt_dist    = [rt_dist; (12-out.data.d11Dist0(out.rt.rt_sel)).*out.data.sgnResp(out.rt.rt_sel)];
        
        rt_res      = [rt_res; prt - out.rt.rt_yhat];
        
        rt_pt       = [rt_pt; out.rt.rt_X(:,1).*pt_count];
        
        rt_tr       = [rt_tr; cumsum(out.rt.rt_X(:,1))<500];
        
        
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
        
        acc_targ5     = [acc_targ5; out.data.d5Targ0(out.acc.acc_sel)];
        acc_dist5     = [acc_dist5; 6-out.data.d5Dist0(out.acc.acc_sel)];
        
        acc_tr       = [acc_tr; cumsum(out.acc.acc_X(:,1))<500];

        
        acc_pt      = [acc_pt; out.acc.acc_X(:,1).*pt_count];
        
        pt_count = pt_count+1;
        
    end
    
    
end


cols = {[36, 123, 160]/255, [242, 95, 92]/255};



rt_tbl = array2table([rt, crt, rt_pred, rt_targ, rt_dist, rt_targ5, rt_dist5, rt_res, rt_pt ], 'VariableNames', ...
    {'rt', 'crt', 'rtY', 'targ', 'dist', 'targ5', 'dist5', 'res', 'pt'});



acc_tbl = array2table([acc, cacc, acc_pred, acc_targ, acc_dist, acc_targ5, acc_dist5, acc_pt], 'VariableNames', ...
    {'acc', 'cacc', 'accY', 'targ', 'dist', 'targ5', 'dist5', 'pt'});


% rt_tbl= rt_tbl(rt_tr==1,:);
% acc_tbl= acc_tbl(acc_tr==1,:);



% ========= Betas
npt = [size(col{1}.data,2), size(col{2}.data,2), size(col{3}.data,2)];


x_col_targ_acc = [1.5, 2.25] - .1;
x_col_targ_rt = [1.5, 2.25] + .1;

x_col_dist_acc = [1, 1.75, 2.5] - .1;
x_col_dist_rt = [1, 1.75, 2.5] + .1;

x_mot_targ_acc = [1.5, 2.25] - .1;
x_mot_targ_rt = [1.5, 2.25] + .1;

x_mot_dist_acc = [1, 1.75, 2.5] - .1;
x_mot_dist_rt = [1, 1.75, 2.5] + .1;


% color


M_col_targ_rt = [...
    getBetaMean(col{2}, 'rt_{targ}'),    ...
    ...
    getBetaMean(col{3}, 'rt_{targ}'),    ...
    ];

CI_col_targ_rt = [...
    getBetaCI(col{2}, 'rt_{targ}'),    ...
    ...
    getBetaCI(col{3}, 'rt_{targ}'),    ...
    ];

all_col_targ_rt2 = getBetaAll(col{2}, 'rt_{targ}');
all_col_targ_rt3 = getBetaAll(col{3}, 'rt_{targ}');


M_col_targ_acc = [...
    getBetaMean(col{2}, 'acc_{targ}'),    ...
    ...
    getBetaMean(col{3}, 'acc_{targ}'),    ...
    ];

CI_col_targ_acc = [...
    getBetaCI(col{2}, 'acc_{targ}'),    ...
    ...
    getBetaCI(col{3}, 'acc_{targ}'),    ...
    ];

all_col_targ_acc2 = getBetaAll(col{2}, 'acc_{targ}');
all_col_targ_acc3 = getBetaAll(col{3}, 'acc_{targ}');



M_col_dist_rt = [...
    getBetaMean(col{1}, 'rt_{dist}'),    ...
    ...
    getBetaMean(col{2}, 'rt_{dist}'),    ...
    ...
    getBetaMean(col{3}, 'rt_{dist}'),    ...
    ];

CI_col_dist_rt = [...
    getBetaCI(col{1}, 'rt_{dist}'),    ...
    ...
    getBetaCI(col{2}, 'rt_{dist}'),    ...
    ...
    getBetaCI(col{3}, 'rt_{dist}'),    ...
    ];

all_col_dist_rt1 = getBetaAll(col{1}, 'rt_{dist}');
all_col_dist_rt2 = getBetaAll(col{2}, 'rt_{dist}');
all_col_dist_rt3 = getBetaAll(col{3}, 'rt_{dist}');



M_col_dist_acc = [...
    getBetaMean(col{1}, 'acc_{dist}'),    ...
    ...
    getBetaMean(col{2}, 'acc_{dist}'),    ...
    ...
    getBetaMean(col{3}, 'acc_{dist}'),    ...
    ];


CI_col_dist_acc = [...
    getBetaCI(col{1}, 'acc_{dist}'),    ...
    ...
    getBetaCI(col{2}, 'acc_{dist}'),    ...
    ...
    getBetaCI(col{3}, 'acc_{dist}'),    ...
    ];

all_col_dist_acc1 = getBetaAll(col{1}, 'acc_{dist}');
all_col_dist_acc2 = getBetaAll(col{2}, 'acc_{dist}');
all_col_dist_acc3 = getBetaAll(col{3}, 'acc_{dist}');





% ============ motion

M_mot_targ_rt = [...
    getBetaMean(mot{2}, 'rt_{targ}'),    ...
    ...
    getBetaMean(mot{3}, 'rt_{targ}'),    ...
    ];

CI_mot_targ_rt = [...
    getBetaCI(mot{2}, 'rt_{targ}'),    ...
    ...
    getBetaCI(mot{3}, 'rt_{targ}'),    ...
    ];

all_mot_targ_rt2 = getBetaAll(mot{2}, 'rt_{targ}');
all_mot_targ_rt3 = getBetaAll(mot{3}, 'rt_{targ}');

M_mot_targ_acc = [...
    getBetaMean(mot{2}, 'acc_{targ}'),    ...
    ...
    getBetaMean(mot{3}, 'acc_{targ}'),    ...
    ];



CI_mot_targ_acc = [...
    getBetaCI(mot{2}, 'acc_{targ}'),    ...
    ...
    getBetaCI(mot{3}, 'acc_{targ}'),    ...
    ];

all_mot_targ_acc2 = getBetaAll(mot{2}, 'acc_{targ}');
all_mot_targ_acc3 = getBetaAll(mot{3}, 'acc_{targ}');


M_mot_dist_rt = [...
    getBetaMean(mot{1}, 'rt_{dist}'),    ...
    ...
    getBetaMean(mot{2}, 'rt_{dist}'),    ...
    ...
    getBetaMean(mot{3}, 'rt_{dist}'),    ...
    ];

CI_mot_dist_rt = [...
    getBetaCI(mot{1}, 'rt_{dist}'),    ...
    ...
    getBetaCI(mot{2}, 'rt_{dist}'),    ...
    ...
    getBetaCI(mot{3}, 'rt_{dist}'),    ...
    ];

all_mot_dist_rt1 = getBetaAll(mot{1}, 'rt_{dist}');
all_mot_dist_rt2 = getBetaAll(mot{2}, 'rt_{dist}');
all_mot_dist_rt3 = getBetaAll(mot{3}, 'rt_{dist}');


M_mot_dist_acc = [...
    getBetaMean(mot{1}, 'acc_{dist}'),    ...
    ...
    getBetaMean(mot{2}, 'acc_{dist}'),    ...
    ...
    getBetaMean(mot{3}, 'acc_{dist}'),    ...
    ];


CI_mot_dist_acc = [...
    getBetaCI(mot{1}, 'acc_{dist}'),    ...
    ...
    getBetaCI(mot{2}, 'acc_{dist}'),    ...
    ...
    getBetaCI(mot{3}, 'acc_{dist}'),    ...
    ];

all_mot_dist_acc1 = getBetaAll(mot{1}, 'acc_{dist}');
all_mot_dist_acc2 = getBetaAll(mot{2}, 'acc_{dist}');
all_mot_dist_acc3 = getBetaAll(mot{3}, 'acc_{dist}');







% plot
figure('units','inch','position',[0,0,2*7.5,2*5.5]); hold on;
tiledlayout(2,20, 'TileSpacing', 'compact', 'Padding', 'compact');





% =================== TARGET


nexttile([1,6]); hold on;

% ============== targ acc ==============
grpTbl = grpstats(acc_tbl, {'pt','targ'});
grpTbl = grpstats(grpTbl, {'targ'}, {'mean','sem'});

yyaxis('left');

% plot([1:length(grpTbl.mean_mean_acc)],(grpTbl.mean_mean_acc),'ob', 'LineWidth', .5);
ep = errorbar(1.15:11.15, grpTbl.mean_mean_acc, grpTbl.sem_mean_cacc,...
    'o', 'color', cols{1}, 'MarkerFaceColor', cols{1}, 'LineWidth', 1);

plot(1.15:11.15, grpTbl.mean_mean_accY, '-', 'color', cols{1}, 'LineWidth', 2);
set(gca, 'TickDir', 'out', 'LineWidth', 1)
yticks([floor(min(ep.YData)*100)/100, ceil(max(ep.YData)*100)/100])

targAcc_ylim = ylim
[mean(targAcc_ylim), range(targAcc_ylim)/2]
ylim([mean(targAcc_ylim)-.07, mean(targAcc_ylim)+.07])
yticks(ylim)

% ============== targ rt ==============
grpTbl = grpstats(rt_tbl, {'pt','targ'}, {'mean'});
grpTbl = grpstats(grpTbl, {'targ'}, {'mean','sem'});

yyaxis('right');

% plot(exp(grpTbl.mean_mean_rt) + mean(all_ndt), 'ob', 'LineWidth', .5);
ep=errorbar(.85:10.85, exp(grpTbl.mean_mean_rt) + mean(all_ndt), exp(grpTbl.mean_mean_rt).*(grpTbl.sem_mean_crt),...
    'o', 'color', cols{2}, 'MarkerFaceColor', cols{2}, 'LineWidth', 1);

plot(.85:10.85, exp(grpTbl.mean_mean_rtY) + mean(all_ndt), '-', 'color', cols{2}, 'LineWidth', 2);

set(gca, 'TickDir', 'out', 'LineWidth', 1)
yticks([floor(min(ep.YData)*100)/100, ceil(max(ep.YData)*100)/100])


title('Target (Attend-Color)')

xlim([.5 11.5])
xticks([])
xticklabels({'low coh', 'high coh'})

xlabel('target coherence')



targRT_ylim = ylim
[mean(targRT_ylim), range(targRT_ylim)/2]
ylim([mean(targRT_ylim)-.06, mean(targRT_ylim)+.06])
yticks(ylim)






% plot
dotOff = 0;

% col targ
n=nexttile([1,4]); hold on;


yyaxis('left');
plot(repmat(x_col_targ_acc(1), [npt(2),1])+randn([npt(2),1]).*.02.*ksdensity(all_col_targ_acc2,all_col_targ_acc2) - dotOff, all_col_targ_acc2,   '.', 'color', cols{1});
plot(repmat(x_col_targ_acc(2), [npt(3),1])+randn([npt(3),1]).*.02.*ksdensity(all_col_targ_acc3,all_col_targ_acc3) - dotOff, all_col_targ_acc3,   '.', 'color', cols{1});
errorbar(x_col_targ_acc+dotOff, M_col_targ_acc,CI_col_targ_acc, 'ok', 'MarkerFaceColor', 'k', 'LineWidth', 2, 'MarkerSize', 7);
set(gca, 'TickDir', 'out', 'LineWidth', 1)
n.YLim = [-6, 6];
yticks([-6, 6])

yyaxis('right');
plot(repmat(x_col_targ_rt(1), [npt(2),1])+randn([npt(2),1]).*.02.*ksdensity(all_col_targ_rt2,all_col_targ_rt2) - dotOff, all_col_targ_rt2,   '.', 'color', cols{2});
plot(repmat(x_col_targ_rt(2), [npt(3),1])+randn([npt(3),1]).*.02.*ksdensity(all_col_targ_rt3,all_col_targ_rt3) - dotOff, all_col_targ_rt3,   '.', 'color', cols{2});
errorbar(x_col_targ_rt+dotOff, M_col_targ_rt, CI_col_targ_rt,  'ok', 'MarkerFaceColor', 'k', 'LineWidth', 2, 'MarkerSize', 7);
set(gca, 'TickDir', 'out', 'LineWidth', 1)
n.YLim = [-1.5, 1.5];
yticks([-1.5, 1.5])

xlim([1,3])
xlim([x_col_targ_acc(1)-.15, x_col_targ_rt(end) + .15])
yline(0, '-k', 'LineWidth', 1);
title('attend-color')

xticks(x_col_targ_acc+.1)
xticklabels({'Exp2', 'Exp3'});






% === motion



mdl = mot;


[rt, crt, rt_pred, rt_targ, rt_dist, rt_res, rt_pt,...
    acc, cacc, acc_pred, acc_targ, acc_dist, acc_res, acc_pt,...
    all_ndt,...
    ] = deal([]);

pt_count = 1;

for mm = 1:length(mot)
    
    
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
    {'rt', 'crt', 'rtY', 'targ', 'dist', 'res', 'pt'});

acc_tbl = array2table([acc, cacc, acc_pred, acc_targ, acc_dist, acc_pt], 'VariableNames', ...
    {'acc', 'cacc', 'accY', 'targ', 'dist', 'pt'});






nexttile([1,6]); hold on;

% ============== targ acc ==============
grpTbl = grpstats(acc_tbl, {'pt','targ'});
grpTbl = grpstats(grpTbl, {'targ'}, {'mean','sem'});

yyaxis('left');

% plot([1:length(grpTbl.mean_mean_acc)],(grpTbl.mean_mean_acc),'ob', 'LineWidth', .5);
ep = errorbar(1.15:11.15, grpTbl.mean_mean_acc, grpTbl.sem_mean_cacc,...
    'o', 'color', cols{1}, 'MarkerFaceColor', cols{1}, 'LineWidth', 1);

plot(1.15:11.15, grpTbl.mean_mean_accY, '-', 'color', cols{1}, 'LineWidth', 2);
set(gca, 'TickDir', 'out', 'LineWidth', 1)
yticks([floor(min(ep.YData)*100)/100, ceil(max(ep.YData)*100)/100])

targAccMot_ylim = ylim
[mean(targAccMot_ylim), range(targAccMot_ylim)/2]

ylim([mean(targAccMot_ylim)-.07, mean(targAccMot_ylim)+.07])
yticks(ylim)

% ============== targ rt ==============
grpTbl = grpstats(rt_tbl, {'pt','targ'}, {'mean'});
grpTbl = grpstats(grpTbl, {'targ'}, {'mean','sem'});

yyaxis('right');

% plot(exp(grpTbl.mean_mean_rt) + mean(all_ndt), 'ob', 'LineWidth', .5);
ep=errorbar(.85:10.85, exp(grpTbl.mean_mean_rt) + mean(all_ndt), exp(grpTbl.mean_mean_rt).*(grpTbl.sem_mean_crt),...
    'o', 'color', cols{2}, 'MarkerFaceColor', cols{2}, 'LineWidth', 1);

plot(.85:10.85, exp(grpTbl.mean_mean_rtY) + mean(all_ndt), '-', 'color', cols{2}, 'LineWidth', 2);

set(gca, 'TickDir', 'out', 'LineWidth', 1)
yticks([floor(min(ep.YData)*100)/100, ceil(max(ep.YData)*100)/100])


title('Target (Attend-Motion)')

xlim([.5 11.5])
xticks([])
xticklabels({'low coh', 'high coh'})

xlabel('target coherence')


targRTMot_ylim = ylim
[mean(targRTMot_ylim), range(targRTMot_ylim)/2]

ylim([mean(targRTMot_ylim)-.06, mean(targRTMot_ylim)+.06])
yticks(ylim)











% mot targ
n=nexttile([1,4]); hold on;


yyaxis('left');
plot(repmat(x_mot_targ_acc(1), [npt(2),1])+randn([npt(2),1]).*.02.*ksdensity(all_mot_targ_acc2,all_mot_targ_acc2) - dotOff, all_mot_targ_acc2,   '.', 'color', cols{1});
plot(repmat(x_mot_targ_acc(2), [npt(3),1])+randn([npt(3),1]).*.02.*ksdensity(all_mot_targ_acc3,all_mot_targ_acc3) - dotOff, all_mot_targ_acc3,   '.', 'color', cols{1});
errorbar(x_mot_targ_acc+dotOff, M_mot_targ_acc, CI_mot_targ_acc, 'ok', 'MarkerFaceColor', 'k', 'LineWidth', 2, 'MarkerSize', 7);
set(gca, 'TickDir', 'out', 'LineWidth', 1)
n.YLim = [-6, 6];
yticks([-6, 6])


yyaxis('right');
plot(repmat(x_mot_targ_rt(1), [npt(2),1])+randn([npt(2),1]).*.02.*ksdensity(all_mot_targ_rt2,all_mot_targ_rt2) - dotOff, all_mot_targ_rt2,   '.', 'color', cols{2});
plot(repmat(x_mot_targ_rt(2), [npt(3),1])+randn([npt(3),1]).*.02.*ksdensity(all_mot_targ_rt3,all_mot_targ_rt3) - dotOff, all_mot_targ_rt3,   '.', 'color', cols{2});
errorbar(x_mot_targ_rt+dotOff, M_mot_targ_rt, CI_mot_targ_rt, 'ok', 'MarkerFaceColor', 'k', 'LineWidth', 2, 'MarkerSize', 7);
set(gca, 'TickDir', 'out', 'LineWidth', 1)
n.YLim = [-1.5, 1.5];
yticks([-1.5, 1.5])



xlim([x_mot_targ_acc(1)-.15, x_mot_targ_rt(end) + .15])
xticks(x_col_targ_acc+.1)
xticklabels({'Exp2', 'Exp3'});
yline(0, '-k', 'LineWidth', 1);
title('attend-motion')




















% =================== DISTRACTOR





mdl = col;


[rt, crt, rt_pred, rt_targ, rt_dist, rt_targ5, rt_dist5, rt_res, rt_pt, rt_tr,...
    acc, cacc, acc_pred, acc_targ, acc_dist, acc_targ5, acc_dist5, acc_res, acc_pt, acc_tr,...
    all_ndt,...
    ] = deal([]);

pt_count = 1;

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
        [~, out]=eval([results.likfun, '(param(tt,:), data(tt), [], [], 0, [])']);
        
        
        prt         = log(out.rt.rt(out.rt.rt_sel) - ndt(tt));
        prt(out.rt.rt(out.rt.rt_sel) <= ndt(tt)) = nan;
        rt          = [rt; prt];
        crt          = [crt; prt - nanmean(prt)];
        
        
        rt_pred     = [rt_pred; out.rt.rt_yhat];
        
        rt_targ     = [rt_targ; out.data.d11Targ0(out.rt.rt_sel)];
        rt_dist     = [rt_dist; 12-out.data.d11Dist0(out.rt.rt_sel)];
        
        rt_targ5     = [rt_targ5; out.data.d5Targ0(out.rt.rt_sel)];
        rt_dist5     = [rt_dist5; 6-out.data.d5Dist0(out.rt.rt_sel)];
        
        %         rt_targ    = [rt_targ; out.data.d11Targ0(out.rt.rt_sel).*out.data.sgnResp(out.rt.rt_sel)];
        %         rt_dist    = [rt_dist; (12-out.data.d11Dist0(out.rt.rt_sel)).*out.data.sgnResp(out.rt.rt_sel)];
        
        rt_res      = [rt_res; prt - out.rt.rt_yhat];
        
        rt_pt       = [rt_pt; out.rt.rt_X(:,1).*pt_count];
        
        rt_tr       = [rt_tr; cumsum(out.rt.rt_X(:,1))<500];
        
        
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
        
        acc_targ5     = [acc_targ5; out.data.d5Targ0(out.acc.acc_sel)];
        acc_dist5     = [acc_dist5; 6-out.data.d5Dist0(out.acc.acc_sel)];
        
        acc_tr       = [acc_tr; cumsum(out.acc.acc_X(:,1))<500];

        
        acc_pt      = [acc_pt; out.acc.acc_X(:,1).*pt_count];
        
        pt_count = pt_count+1;
        
    end
    
    
end


cols = {[36, 123, 160]/255, [242, 95, 92]/255};



rt_tbl = array2table([rt, crt, rt_pred, rt_targ, rt_dist, rt_targ5, rt_dist5, rt_res, rt_pt ], 'VariableNames', ...
    {'rt', 'crt', 'rtY', 'targ', 'dist', 'targ5', 'dist5', 'res', 'pt'});



acc_tbl = array2table([acc, cacc, acc_pred, acc_targ, acc_dist, acc_targ5, acc_dist5, acc_pt], 'VariableNames', ...
    {'acc', 'cacc', 'accY', 'targ', 'dist', 'targ5', 'dist5', 'pt'});


% rt_tbl= rt_tbl(rt_tr==1,:);
% acc_tbl= acc_tbl(acc_tr==1,:);











nexttile([1,6]); hold on;

% ============== dist acc ==============
yyaxis('left');

grpTbl = grpstats(acc_tbl, {'pt','dist'});
grpTbl = grpstats(grpTbl, {'dist'}, {'mean','sem'});

% plot([1:length(grpTbl.mean_mean_acc)], (grpTbl.mean_mean_acc),'ok', 'LineWidth', .5);
ep=errorbar(1.15:11.15, grpTbl.mean_mean_acc, grpTbl.sem_mean_cacc,...
    'o', 'color', cols{1}, 'MarkerFaceColor', cols{1}, 'LineWidth', 1);
plot(1.15:11.15, grpTbl.mean_mean_accY, '-', 'color', cols{1}, 'LineWidth', 2);
set(gca, 'TickDir', 'out', 'LineWidth', 1)
yticks([floor(min(ep.YData)*100)/100, ceil(max(ep.YData)*100)/100])
xticks([])

distAcc_ylim = ylim
[mean(distAcc_ylim), range(distAcc_ylim)/2]
ylim([mean(distAcc_ylim)-.05, mean(distAcc_ylim)+.05])
yticks(ylim)

% ============== dist rt ==============
yyaxis('right');

grpTbl = grpstats(rt_tbl, {'pt','dist'}, {'mean'});
grpTbl = grpstats(grpTbl, {'dist'}, {'mean','sem'});



% plot(exp(grpTbl.mean_mean_rt) + mean(all_ndt), 'ok', 'LineWidth', .5);
ep = errorbar(.85:10.85, exp(grpTbl.mean_mean_rt) + mean(all_ndt), exp(grpTbl.mean_mean_rt).*(grpTbl.sem_mean_crt),...
    'o', 'color', cols{2}, 'MarkerFaceColor', cols{2}, 'LineWidth', 1);
plot(.85:10.85, exp(grpTbl.mean_mean_rtY)+ mean(all_ndt), '-', 'color', cols{2}, 'LineWidth', 2);
set(gca, 'TickDir', 'out', 'LineWidth', 1)
title('Distractor (Attend-Color)')
% ylim(targRange);

yticks([floor(min(ep.YData)*100)/100, ceil(max(ep.YData)*100)/100])
xlim([.5 11.5])
xticks([])
xlabel('distractor congruence')

distRT_ylim = ylim
[mean(distRT_ylim), range(distRT_ylim)/2]
ylim([mean(distRT_ylim)-.03, mean(distRT_ylim)+.03])
yticks(ylim)






% col dist
n=nexttile([1,4]); hold on;



yyaxis('left');
plot(repmat(x_col_dist_acc(1), [npt(1),1])+randn([npt(1),1]).*.01.*ksdensity(all_col_dist_acc1,all_col_dist_acc1) - dotOff, all_col_dist_acc1,   '.', 'color', cols{1});
plot(repmat(x_col_dist_acc(2), [npt(2),1])+randn([npt(2),1]).*.01.*ksdensity(all_col_dist_acc2,all_col_dist_acc2) - dotOff, all_col_dist_acc2,   '.', 'color', cols{1});
plot(repmat(x_col_dist_acc(3), [npt(3),1])+randn([npt(3),1]).*.01.*ksdensity(all_col_dist_acc3,all_col_dist_acc3) - dotOff, all_col_dist_acc3,   '.', 'color', cols{1});
errorbar(x_col_dist_acc+dotOff, M_col_dist_acc, CI_col_dist_acc, 'ok', 'MarkerFaceColor', 'k', 'LineWidth', 2, 'MarkerSize', 7);
set(gca, 'TickDir', 'out', 'LineWidth', 1)
n.YLim = [-1.8, 1.8];
yticks([-1.8, 1.8])

yyaxis('right');
plot(repmat(x_col_dist_rt(1), [npt(1),1])+randn([npt(1),1]).*.0025.*ksdensity(all_col_dist_rt1,all_col_dist_rt1) - dotOff, all_col_dist_rt1,   '.', 'color', cols{2});
plot(repmat(x_col_dist_rt(2), [npt(2),1])+randn([npt(2),1]).*.0025.*ksdensity(all_col_dist_rt2,all_col_dist_rt2) - dotOff, all_col_dist_rt2,   '.', 'color', cols{2});
plot(repmat(x_col_dist_rt(3), [npt(3),1])+randn([npt(3),1]).*.0025.*ksdensity(all_col_dist_rt3,all_col_dist_rt3) - dotOff, all_col_dist_rt3,   '.', 'color', cols{2});
errorbar(x_col_dist_rt+dotOff, M_col_dist_rt, CI_col_dist_rt, 'ok', 'MarkerFaceColor', 'k', 'LineWidth', 2, 'MarkerSize', 7);
set(gca, 'TickDir', 'out', 'LineWidth', 1)
n.YLim = [-.15, .15];
yticks([-.15, .15])




xlim([x_col_dist_acc(1)-.15, x_col_dist_rt(end) + .15])
xticks(x_col_dist_acc+.1)
xticklabels({'Exp1', 'Exp2', 'Exp3'});

yline(0, '-k', 'LineWidth', 1);

title('attend-color')









% === motion



mdl = mot;


[rt, crt, rt_pred, rt_targ, rt_dist, rt_res, rt_pt,...
    acc, cacc, acc_pred, acc_targ, acc_dist, acc_res, acc_pt,...
    all_ndt,...
    ] = deal([]);

pt_count = 1;

for mm = 1:length(mot)
    
    
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
    {'rt', 'crt', 'rtY', 'targ', 'dist', 'res', 'pt'});

acc_tbl = array2table([acc, cacc, acc_pred, acc_targ, acc_dist, acc_pt], 'VariableNames', ...
    {'acc', 'cacc', 'accY', 'targ', 'dist', 'pt'});








nexttile([1,6]); hold on;

% ============== dist acc ==============
yyaxis('left');

grpTbl = grpstats(acc_tbl, {'pt','dist'});
grpTbl = grpstats(grpTbl, {'dist'}, {'mean','sem'});

% plot([1:length(grpTbl.mean_mean_acc)], (grpTbl.mean_mean_acc),'ok', 'LineWidth', .5);
ep=errorbar(1.15:11.15, grpTbl.mean_mean_acc, grpTbl.sem_mean_cacc,...
    'o', 'color', cols{1}, 'MarkerFaceColor', cols{1}, 'LineWidth', 1);
plot(1.15:11.15, grpTbl.mean_mean_accY, '-', 'color', cols{1}, 'LineWidth', 2);
set(gca, 'TickDir', 'out', 'LineWidth', 1)
yticks([floor(min(ep.YData)*100)/100, ceil(max(ep.YData)*100)/100])
xticks([])

distAccMot_ylim = ylim
[mean(distAccMot_ylim), range(distAccMot_ylim)/2]
ylim([mean(distAccMot_ylim)-.05, mean(distAccMot_ylim)+.05])
yticks(ylim)

% ============== dist rt ==============
yyaxis('right');

grpTbl = grpstats(rt_tbl, {'pt','dist'}, {'mean'});
grpTbl = grpstats(grpTbl, {'dist'}, {'mean','sem'});



% plot(exp(grpTbl.mean_mean_rt) + mean(all_ndt), 'ok', 'LineWidth', .5);
ep = errorbar(.85:10.85, exp(grpTbl.mean_mean_rt) + mean(all_ndt), exp(grpTbl.mean_mean_rt).*(grpTbl.sem_mean_crt),...
    'o', 'color', cols{2}, 'MarkerFaceColor', cols{2}, 'LineWidth', 1);
plot(.85:10.85, exp(grpTbl.mean_mean_rtY)+ mean(all_ndt), '-', 'color', cols{2}, 'LineWidth', 2);
set(gca, 'TickDir', 'out', 'LineWidth', 1)
title('Distractor (Attend-Motion)')
% ylim(targRange);

yticks([floor(min(ep.YData)*100)/100, ceil(max(ep.YData)*100)/100])

xlim([.5 11.5])
xticks([])
xlabel('distractor congruence')

distRTMot_ylim = ylim
[mean(distRTMot_ylim), range(distRTMot_ylim)/2]
ylim([mean(distRTMot_ylim)-.03, mean(distRTMot_ylim)+.03])
yticks(ylim)








% mot dist
n=nexttile([1,4]); hold on;


yyaxis('left');
plot(repmat(x_mot_dist_acc(1), [npt(1),1])+randn([npt(1),1]).*.01.*ksdensity(all_mot_dist_acc1,all_mot_dist_acc1) - dotOff, all_mot_dist_acc1,   '.', 'color', cols{1});
plot(repmat(x_mot_dist_acc(2), [npt(2),1])+randn([npt(2),1]).*.01.*ksdensity(all_mot_dist_acc2,all_mot_dist_acc2) - dotOff, all_mot_dist_acc2,   '.', 'color', cols{1});
plot(repmat(x_mot_dist_acc(3), [npt(3),1])+randn([npt(3),1]).*.01.*ksdensity(all_mot_dist_acc3,all_mot_dist_acc3) - dotOff, all_mot_dist_acc3,   '.', 'color', cols{1});
errorbar(x_mot_dist_acc + dotOff, M_mot_dist_acc, CI_mot_dist_acc, 'ok', 'MarkerFaceColor', 'k', 'LineWidth', 2, 'MarkerSize', 7);
set(gca, 'TickDir', 'out', 'LineWidth', 1)
n.YLim = [-1.8, 1.8];
yticks([-1.8, 1.8])


yyaxis('right');
plot(repmat(x_mot_dist_rt(1), [npt(1),1])+randn([npt(1),1]).*.001.*ksdensity(all_mot_dist_rt1,all_mot_dist_rt1) - dotOff, all_mot_dist_rt1,   '.', 'color', cols{2});
plot(repmat(x_mot_dist_rt(2), [npt(2),1])+randn([npt(2),1]).*.001.*ksdensity(all_mot_dist_rt2,all_mot_dist_rt2) - dotOff, all_mot_dist_rt2,   '.', 'color', cols{2});
plot(repmat(x_mot_dist_rt(3), [npt(3),1])+randn([npt(3),1]).*.001.*ksdensity(all_mot_dist_rt3,all_mot_dist_rt3) - dotOff, all_mot_dist_rt3,   '.', 'color', cols{2});
errorbar(x_mot_dist_rt + dotOff, M_mot_dist_rt, CI_mot_dist_rt, 'ok', 'MarkerFaceColor', 'k', 'LineWidth', 2, 'MarkerSize', 7);
set(gca, 'TickDir', 'out', 'LineWidth', 1)

n.YLim = [-.15, .15];
yticks([-.15, .15])




xlim([x_mot_dist_acc(1)-.15, x_mot_dist_rt(end) + .15])
xticks(x_mot_dist_acc+.1)
xticklabels({'Exp1', 'Exp2', 'Exp3'});

yline(0, '-k', 'LineWidth', 1);

title('attend-motion')
















% Plot interaction ========================================================





mdl = col;


[rt, crt, rt_pred, rt_targ, rt_dist, rt_targ5, rt_dist5, rt_res, rt_pt,...
    acc, cacc, acc_pred, acc_targ, acc_dist, acc_targ5, acc_dist5, acc_res, acc_pt,...
    all_ndt,...
    ] = deal([]);

pt_count = 1;

for mm = 2:3
    
    
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
        
        rt_targ5     = [rt_targ5; out.data.d5Targ0(out.rt.rt_sel)];
        rt_dist5     = [rt_dist5; 6-out.data.d5Dist0(out.rt.rt_sel)];
        
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
        
        acc_targ5     = [acc_targ5; out.data.d5Targ0(out.acc.acc_sel)];
        acc_dist5     = [acc_dist5; 6-out.data.d5Dist0(out.acc.acc_sel)];
        
        
        acc_pt      = [acc_pt; out.acc.acc_X(:,1).*pt_count];
        
        pt_count = pt_count+1;
        
    end
    
    
end


cols = {[36, 123, 160]/255, [242, 95, 92]/255};



rt_tbl = array2table([rt, crt, rt_pred, rt_targ, rt_dist, rt_targ5, rt_dist5, rt_res, rt_pt ], 'VariableNames', ...
    {'rt', 'crt', 'rtY', 'targ', 'dist', 'targ5', 'dist5', 'res', 'pt'});

acc_tbl = array2table([acc, cacc, acc_pred, acc_targ, acc_dist, acc_targ5, acc_dist5, acc_pt], 'VariableNames', ...
    {'acc', 'cacc', 'accY', 'targ', 'dist', 'targ5', 'dist5', 'pt'});




mdl = noX;


[rt, crt, rt_pred, rt_targ, rt_dist, rt_targ5, rt_dist5, rt_res, rt_pt,...
    acc, cacc, acc_pred, acc_targ, acc_dist, acc_targ5, acc_dist5, acc_res, acc_pt,...
    all_ndtx,...
    ] = deal([]);

pt_count = 1;

for mm = 2:3
    
    
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
        
        rt_targ5     = [rt_targ5; out.data.d5Targ0(out.rt.rt_sel)];
        rt_dist5     = [rt_dist5; 6-out.data.d5Dist0(out.rt.rt_sel)];
        
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
        
        acc_targ5     = [acc_targ5; out.data.d5Targ0(out.acc.acc_sel)];
        acc_dist5     = [acc_dist5; 6-out.data.d5Dist0(out.acc.acc_sel)];
        
        
        acc_pt      = [acc_pt; out.acc.acc_X(:,1).*pt_count];
        
        pt_count = pt_count+1;
        
    end
    
    
end


cols = {[36, 123, 160]/255, [242, 95, 92]/255};


rt_tbl.rt = rt;
acc_tbl.acc = acc;

rt_tbl.rtYx = rt_pred;
acc_tbl.accYx = acc_pred;











figure('units','inch','position',[0,0,9,9]); hold on;
tiledlayout(2,2, 'TileSpacing', 'compact', 'Padding', 'compact');



parcols = colormap('parula');
colidx = round(linspace(25,225,5));
xoff = linspace(-.2, .2, 5);






% ============ distractor sorted ============ 
% Accuracy
nexttile; hold on;
grpTbl = grpstats(acc_tbl, {'pt','targ5', 'dist5'});
grpTbl = grpstats(grpTbl, {'targ5', 'dist5'}, {'mean','sem'});

for ii = 1:5
ep = errorbar([1:5] + xoff(ii), grpTbl.mean_mean_acc(grpTbl.dist5 == ii), grpTbl.sem_mean_cacc(grpTbl.dist5 == ii),...
    'o', 'color', parcols(colidx(ii),:), 'MarkerFaceColor', parcols(colidx(ii),:), 'LineWidth', 1);

plot([1:5] + xoff(ii), grpTbl.mean_mean_accY(grpTbl.dist5 == ii), '-', 'color', parcols(colidx(ii),:), 'LineWidth', 2);
plot([1:5] + xoff(ii), grpTbl.mean_mean_accYx(grpTbl.dist5 == ii), '--', 'color', parcols(colidx(ii),:), 'LineWidth', 2);

end


set(gca, 'TickDir', 'out', 'LineWidth', 1)
yticks(ylim)
xlim([.5 5.5])
xticks([])
xlabel('Target Coherence')
title('Distractor-Sorted Accuracy')
ylabel('Accuracy')



% RT
grpTbl = grpstats(rt_tbl, {'pt','targ5', 'dist5'}, {'mean'});
grpTbl = grpstats(grpTbl, {'targ5', 'dist5'}, {'mean','sem'});

nexttile; hold on;

for ii = 1:5
ep=errorbar([1:5] + xoff(ii), exp(grpTbl.mean_mean_rt(grpTbl.dist5 == ii)) + mean([all_ndt]), exp(grpTbl.mean_mean_rt(grpTbl.dist5 == ii)).*(grpTbl.sem_mean_crt(grpTbl.dist5 == ii)),...
    'o', 'color', parcols(colidx(ii),:), 'MarkerFaceColor', parcols(colidx(ii),:), 'LineWidth', 1);

plot([1:5] + xoff(ii), exp(grpTbl.mean_mean_rtY(grpTbl.dist5 == ii)) + mean([all_ndt]), '-', 'color', parcols(colidx(ii),:), 'LineWidth', 2);
plot([1:5] + xoff(ii), exp(grpTbl.mean_mean_rtYx(grpTbl.dist5 == ii)) + mean([all_ndt]), '--', 'color', parcols(colidx(ii),:), 'LineWidth', 2);

end

set(gca, 'TickDir', 'out', 'LineWidth', 1)
yticks(ylim)
xlim([.5 5.5])
xticks([])

xlabel('Target Coherence')
title('Distractor-Sorted RT')
ylabel('Reaction Time')




% ============ target sorted ============ 
% Accuracy
nexttile; hold on;
grpTbl = grpstats(acc_tbl, {'pt','targ5', 'dist5'});
grpTbl = grpstats(grpTbl, {'targ5', 'dist5'}, {'mean','sem'});


for ii = 1:5
ep = errorbar([1:5] + xoff(ii), grpTbl.mean_mean_acc(grpTbl.targ5 == ii), grpTbl.sem_mean_cacc(grpTbl.targ5 == ii),...
    'o', 'color', parcols(colidx(ii),:), 'MarkerFaceColor', parcols(colidx(ii),:), 'LineWidth', 1);

plot([1:5] + xoff(ii), grpTbl.mean_mean_accY(grpTbl.targ5 == ii), '-', 'color', parcols(colidx(ii),:), 'LineWidth', 2);
plot([1:5] + xoff(ii), grpTbl.mean_mean_accYx(grpTbl.targ5 == ii), '--', 'color', parcols(colidx(ii),:), 'LineWidth', 2);

end


set(gca, 'TickDir', 'out', 'LineWidth', 1)
yticks(ylim)
xlim([.5 5.5])
xticks([])
xlabel('Distractor Congruence')
title('Target-Sorted Accuracy')
ylabel('Accuracy')



% RT
grpTbl = grpstats(rt_tbl, {'pt','targ5', 'dist5'}, {'mean'});
grpTbl = grpstats(grpTbl, {'targ5', 'dist5'}, {'mean','sem'});

nexttile; hold on;

for ii = 1:5
ep=errorbar([1:5] + xoff(ii), exp(grpTbl.mean_mean_rt(grpTbl.targ5 == ii)) + mean([all_ndt]), exp(grpTbl.mean_mean_rt(grpTbl.targ5 == ii)).*(grpTbl.sem_mean_crt(grpTbl.targ5 == ii)),...
    'o', 'color', parcols(colidx(ii),:), 'MarkerFaceColor', parcols(colidx(ii),:), 'LineWidth', 1);

plot([1:5] + xoff(ii), exp(grpTbl.mean_mean_rtY(grpTbl.targ5 == ii)) + mean([all_ndt]), '-', 'color', parcols(colidx(ii),:), 'LineWidth', 2);
plot([1:5] + xoff(ii), exp(grpTbl.mean_mean_rtYx(grpTbl.targ5 == ii)) + mean([all_ndt]), '--', 'color', parcols(colidx(ii),:), 'LineWidth', 2);

end

set(gca, 'TickDir', 'out', 'LineWidth', 1)
yticks(ylim)
xlim([.5 5.5])
xticks([])

xlabel('Distractor Congruence')
ylabel('Reaction Time')
title('Target-Sorted RT')


%% targ-dist interaction model comparision
figure; hold on;

nexttile; hold on;
violinplot([noX{2}.results.bf.bic - col{2}.results.bf.bic]);
ylim([-25, 25]);
yline([0])
title('exp2 BIC')

[~,~,~,pxp,bor] = spm_BMS(-[noX{2}.results.bf.bic; col{2}.results.bf.bic]')

nexttile; hold on;
violinplot([noX{2}.results.bf.aic - col{2}.results.bf.aic]);
ylim([-25, 25]);
yline([0])
title('exp2 AIC')

[~,~,~,pxp,bor] = spm_BMS(-[noX{2}.results.bf.aic; col{2}.results.bf.aic]')


nexttile; hold on;
violinplot([noX{3}.results.bf.bic - col{3}.results.bf.bic]);
ylim([-25, 25]);
yline([0])
title('exp3 BIC')

[~,~,~,pxp,bor] = spm_BMS(-[noX{3}.results.bf.bic; col{3}.results.bf.bic]')

nexttile; hold on;
violinplot([noX{3}.results.bf.aic - col{3}.results.bf.aic]);
ylim([-25, 25]);
yline([0])
title('exp3 AIC')

[~,~,~,pxp,bor] = spm_BMS(-[noX{3}.results.bf.aic; col{3}.results.bf.aic]')






%% ====== DYNAMICS =========================================================
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


%% plot dynamic



load('batlow.mat');
cols = batlow;
cIdx = round(linspace(1, length(cols), 7));


rtacc_cols = {[36, 123, 160]/255, [242, 95, 92]/255};


clear mdl
mdl = col;


[rt, crt, d5rt, rt_pred, rt_targ, rt_dist, rt_acc, rt_pt,...
    acc, cacc, acc_pred, acc_targ, acc_dist, acc_rt, acc_pt,...
    all_ndt, acc_pdf,...
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
        crt          = [crt; prt - nanmean(prt)];
        
        
        rt_pred     = [rt_pred; out.rt.rt_yhat];
        
        rt_targ     = [acc_targ; discretize(tiedrank(out.data.targ0(out.rt.rt_sel)),nCohDisc)];
        rt_dist     = [acc_dist; discretize(tiedrank(out.data.dist0(out.rt.rt_sel)),nCohDisc)];
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
        
        acc_rt          = [acc_rt; out.data.rt0(out.acc.acc_sel)];
        acc_pt      = [acc_pt; out.acc.acc_X(:,1).*pt_count];
        
        acc_pdf     = [acc_pdf; out.acc.acc_pdf];
        
        pt_count = pt_count+1;
        
    end
    
    
end



rt_tbl = array2table([rt, crt, rt_pred, rt_targ, rt_dist, rt_acc, rt_pt ], 'VariableNames', ...
    {'rt', 'crt', 'rtY', 'targ', 'dist', 'acc', 'pt'});

acc_tbl = array2table([acc, cacc, acc_pred, acc_rt, d5rt, acc_targ, acc_dist, acc_pdf, acc_pt], 'VariableNames', ...
    {'acc', 'cacc', 'accY', 'rt', 'd5rt', 'targ', 'dist', 'pdf', 'pt'});






% acc_tbl = array2table([acc, acc_pred, rt, d5rt, acc_targ, acc_dist, acc_pt], 'VariableNames', ...
%     {'acc', 'accY','rt', 'd5rt', 'targ', 'dist' 'pt'});


npt = [size(col{1}.data,2), size(col{2}.data,2), size(col{3}.data,2)];
nptSel = [ones(npt(1),1)*1; ones(npt(2),1)*2; ones(npt(3),1)*3];








% ===== start plot
figure('units','inch','position',[0,0,2*5.5,2*5.5]); hold on;
tiledlayout(3,20, 'TileSpacing', 'compact', 'Padding', 'compact');








% ============== plot targ dynamics acc ==============
grpTbl = grpstats(acc_tbl, {'pt', 'd5rt', 'targ'});
grpTbl = grpstats(grpTbl, {'d5rt', 'targ'}, {'mean','sem'});


nexttile([1,8]); hold on;

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

xticks([]);


yticks([ceil(min(allY)*100)/100, floor(max(allY)*100)/100])
yticks(ylim)

title('Target \times RT')
xlabel('Target Coherence')

ylabel('relative accuracy')







% ============== targ dynamics RT ==============
grpTbl = grpstats(rt_tbl, {'pt', 'acc', 'targ'});
grpTbl = grpstats(grpTbl, {'acc', 'targ'}, {'mean','sem'});

rt_cIdx = [180, 100];

xoff = [-.15, .15];

n=nexttile([1,8]); hold on;
allY = [];
for ii = 0:1
    
    
    e=errorbar([1:nCohDisc]  + xoff(ii+1), (exp(grpTbl.mean_mean_rt(grpTbl.acc == ii)) + mean(all_ndt)), exp(grpTbl.mean_mean_rt(grpTbl.acc == ii)).*grpTbl.sem_mean_crt(grpTbl.acc == ii),...
        'o', 'LineWidth', 1);
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
yticks(ylim)

ylabel('RT')









% ========================== BETAS

n=nexttile([1,4]); hold on;


% ============== targ dynamics acc betas ==============
x_col_targ_acc = [1.25, 2];


M_col_targ_acc = [...
    getBetaMean(col{2}, 'acc_{targCrt}'),    ...
    ...
    getBetaMean(col{3}, 'acc_{targCrt}'),    ...
    ];

CI_col_targ_acc = [...
    getBetaCI(col{2}, 'acc_{targCrt}'),    ...
    ...
    getBetaCI(col{3}, 'acc_{targCrt}'),    ...
    ];

all_col_targ_acc2 = getBetaAll(col{2}, 'acc_{targCrt}');
all_col_targ_acc3 = getBetaAll(col{3}, 'acc_{targCrt}');


% n=nexttile([1,4]); hold on;
yyaxis('left');


plot(repmat(x_col_targ_acc(1), [npt(2),1])+randn([npt(2),1]).*.4.*ksdensity(all_col_targ_acc2,all_col_targ_acc2), all_col_targ_acc2, '.', 'color', rtacc_cols{1});
plot(repmat(x_col_targ_acc(2), [npt(3),1])+randn([npt(3),1]).*.4.*ksdensity(all_col_targ_acc3,all_col_targ_acc3), all_col_targ_acc3, '.', 'color', rtacc_cols{1});
errorbar(x_col_targ_acc, M_col_targ_acc, CI_col_targ_acc,  'ok', 'MarkerFaceColor', 'k', 'LineWidth', 2, 'MarkerSize', 7);


% n.YLim = [0, Inf];
set(gca, 'TickDir', 'out', 'LineWidth', 1)
% xlim([x_col_targ_acc(1)-.15, x_col_targ_acc(end) + .15])
yline(0, '-k', 'LineWidth', 1);
ylim([-10, 15]);
yticks([-10, 0, 15]);






% ============== targ dynamics RT betas ==============
x_col_targ_acc = [1.5, 2.25];


M_col_targ_acc = [...
    getBetaMean(col{2}, 'rt_{targAcc}'),    ...
    ...
    getBetaMean(col{3}, 'rt_{targAcc}'),    ...
    ];

CI_col_targ_acc = [...
    getBetaCI(col{2}, 'rt_{targAcc}'),    ...
    ...
    getBetaCI(col{3}, 'rt_{targAcc}'),    ...
    ];

all_col_targ_acc2 = getBetaAll(col{2}, 'rt_{targAcc}');
all_col_targ_acc3 = getBetaAll(col{3}, 'rt_{targAcc}');


yyaxis('right');

plot(repmat(x_col_targ_acc(1), [npt(2),1])+randn([npt(2),1]).*.015.*ksdensity(all_col_targ_acc2,all_col_targ_acc2), all_col_targ_acc2, '.', 'color', rtacc_cols{2});
plot(repmat(x_col_targ_acc(2), [npt(3),1])+randn([npt(3),1]).*.015.*ksdensity(all_col_targ_acc3,all_col_targ_acc3), all_col_targ_acc3, '.', 'color', rtacc_cols{2});
errorbar(x_col_targ_acc, M_col_targ_acc, CI_col_targ_acc,  'ok', 'MarkerFaceColor', 'k', 'LineWidth', 2, 'MarkerSize', 7);


ylim([-(.85/1.5), .85]);
yticks([-.56, 0, .85]);

set(gca, 'TickDir', 'out', 'LineWidth', 1)
xlim([1.25-.15, 2.25 + .15])
yline(0, '-k', 'LineWidth', 1);

xticks([1.375, 2.125])

xticklabels({'Exp2', 'Exp 3'})


title('Target Betas')













% ============ Distractor









% ============== dist dynamics acc ==============
grpTbl = grpstats(acc_tbl, {'pt', 'd5rt', 'dist'});
grpTbl = grpstats(grpTbl, {'d5rt', 'dist'}, {'mean','sem'});

nexttile([1,8]); hold on;

xoff = linspace(-.2,.2,5);

yline(0, '-k', 'LineWidth', 1);


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

yticks([ceil(min(allY)*100)/100, floor(max(allY)*100)/100])
yticks(ylim)











% ============== dist dynamics RT ==============
grpTbl = grpstats(rt_tbl, {'pt', 'acc', 'dist'});
grpTbl = grpstats(grpTbl, {'acc', 'dist'}, {'mean','sem'});

xoff = [-.15, .15];


n=nexttile([1,8]); hold on;
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
yticks(ylim)

ylabel('RT')
title('Distractor \times Accuracy')










% ============== dist dynamics acc betas ==============
x_col_targ_acc = [1.5, 2.25];
x_col_dist_acc = [1, 1.75, 2.5]-.125;


M_col_dist_acc = [...
    getBetaMean(col{1}, 'acc_{distCrt}'),    ...
    ...
    getBetaMean(col{2}, 'acc_{distCrt}'),    ...
    ...
    getBetaMean(col{3}, 'acc_{distCrt}'),    ...
    ];

CI_col_dist_acc = [...
    getBetaCI(col{1}, 'acc_{distCrt}'),    ...
    ...
    getBetaCI(col{2}, 'acc_{distCrt}'),    ...
    ...
    getBetaCI(col{3}, 'acc_{distCrt}'),    ...
    ];

all_col_dist_acc1 = getBetaAll(col{1}, 'acc_{distCrt}');
all_col_dist_acc2 = getBetaAll(col{2}, 'acc_{distCrt}');
all_col_dist_acc3 = getBetaAll(col{3}, 'acc_{distCrt}');


n=nexttile([1,4]); hold on;
yyaxis('left');


plot(repmat(x_col_dist_acc(1), [npt(1),1])+randn([npt(1),1]).*.2.*ksdensity(all_col_dist_acc1,all_col_dist_acc1), all_col_dist_acc1,   '.', 'color', rtacc_cols{1});
plot(repmat(x_col_dist_acc(2), [npt(2),1])+randn([npt(2),1]).*.2.*ksdensity(all_col_dist_acc2,all_col_dist_acc2), all_col_dist_acc2,   '.', 'color', rtacc_cols{1});
plot(repmat(x_col_dist_acc(3), [npt(3),1])+randn([npt(3),1]).*.2.*ksdensity(all_col_dist_acc3,all_col_dist_acc3), all_col_dist_acc3,   '.', 'color', rtacc_cols{1});
errorbar(x_col_dist_acc, M_col_dist_acc, CI_col_dist_acc,  'ok', 'MarkerFaceColor', 'k', 'LineWidth', 2, 'MarkerSize', 7);


% n.YLim = [0, Inf];
set(gca, 'TickDir', 'out', 'LineWidth', 1)
xlim([1-.125-.25, 2.5+.125+.25])
yline(0, '-k', 'LineWidth', 1);
ylim([-8, 8]);
yticks([-8,0, 8]);




% ============== dist dynamics RT betas ==============
x_col_targ_acc = [1.5, 2.25];
x_col_dist_acc = [1, 1.75, 2.5]+.125;


M_col_dist_acc = [...
    getBetaMean(col{1}, 'rt_{distAcc}'),    ...
    ...
    getBetaMean(col{2}, 'rt_{distAcc}'),    ...
    ...
    getBetaMean(col{3}, 'rt_{distAcc}'),    ...
    ];

CI_col_dist_acc = [...
    getBetaCI(col{1}, 'rt_{distAcc}'),    ...
    ...
    getBetaCI(col{2}, 'rt_{distAcc}'),    ...
    ...
    getBetaCI(col{3}, 'rt_{distAcc}'),    ...
    ];

all_col_dist_acc1 = getBetaAll(col{1}, 'rt_{distAcc}');
all_col_dist_acc2 = getBetaAll(col{2}, 'rt_{distAcc}');
all_col_dist_acc3 = getBetaAll(col{3}, 'rt_{distAcc}');


% n=nexttile([1,2]); hold on;
yyaxis('right');


plot(repmat(x_col_dist_acc(1), [npt(1),1])+randn([npt(1),1]).*.01.*ksdensity(all_col_dist_acc1,all_col_dist_acc1), all_col_dist_acc1,   '.', 'color', rtacc_cols{2});
plot(repmat(x_col_dist_acc(2), [npt(2),1])+randn([npt(2),1]).*.01.*ksdensity(all_col_dist_acc2,all_col_dist_acc2), all_col_dist_acc2,   '.', 'color', rtacc_cols{2});
plot(repmat(x_col_dist_acc(3), [npt(3),1])+randn([npt(3),1]).*.01.*ksdensity(all_col_dist_acc3,all_col_dist_acc3), all_col_dist_acc3,   '.', 'color', rtacc_cols{2});
errorbar(x_col_dist_acc, M_col_dist_acc, CI_col_dist_acc,  'ok', 'MarkerFaceColor', 'k', 'LineWidth', 2, 'MarkerSize', 7);


% n.YLim = [0, Inf];
set(gca, 'TickDir', 'out', 'LineWidth', 1)
% xlim([x_col_dist_acc(1)-.25, x_col_dist_acc(end) + .25])
yline(0, '-k', 'LineWidth', 1);
ylim([-.5, .5]);
yticks([-.5, 0, .5]);

xticks([1, 1.75, 2.5])
xticklabels({'Exp1', 'Exp2', 'Exp3'})

title('Distractor Betas')















% == plot dynamics LINEAR     ~~~~~~~~~ COLOR ~~~~~~~~~
grpTbl = grpstats(acc_tbl, {'pt'});
% rtMx = repmat([prctile(acc_tbl.rt,5):.001:prctile(acc_tbl.rt,95)], [size(grpTbl,1),1]) - grpTbl.mean_rt;
%
% mRT   = mean(acc_tbl.rt);
% rtRange = (prctile(acc_tbl.rt,5)-mRT):.001:(prctile(acc_tbl.rt,95)-mRT);






rtRange = .5:.001:.9;

clear mRT
for pp = 1:size(grpTbl,1)
    
    mRT(pp,1) = mean(acc_tbl.rt(acc_pt==pp));
    
end
rtMx = repmat(rtRange, [size(grpTbl,1),1]) - mRT;






% estimate
% dyn_targ  = [...
%     getBetaAll(col{1}, 'acc_{targ}') + rtRange.*getBetaAll(col{1}, 'acc_{targCrt}');...
%     getBetaAll(col{2}, 'acc_{targ}') + rtRange.*getBetaAll(col{2}, 'acc_{targCrt}');...
%     getBetaAll(col{3}, 'acc_{targ}') + rtRange.*getBetaAll(col{3}, 'acc_{targCrt}')];
%
% dyn_dist = [...
%     getBetaAll(col{1}, 'acc_{dist}') + rtRange.*getBetaAll(col{1}, 'acc_{distCrt}');...
%     getBetaAll(col{2}, 'acc_{dist}') + rtRange.*getBetaAll(col{2}, 'acc_{distCrt}');...
%     getBetaAll(col{3}, 'acc_{dist}') + rtRange.*getBetaAll(col{3}, 'acc_{distCrt}')];



dyn_lapse  = [...
    1./(1+exp(-(getBetaAll(col{2}, 'acc_{lapse}') + rtMx(nptSel==2,:).*getBetaAll(col{2}, 'acc_{crt}'))));... + (rtMx(nptSel==2,:).^2).*getBetaAll(col{2}, 'acc_{crt2}'))));
    1./(1+exp(-(getBetaAll(col{3}, 'acc_{lapse}') + rtMx(nptSel==3,:).*getBetaAll(col{3}, 'acc_{crt}'))));... + (rtMx(nptSel==3,:).^2).*getBetaAll(col{3}, 'acc_{crt2}'))));
    ];




dyn_targ  = (1-dyn_lapse).*[...
    getBetaAll(col{2}, 'acc_{targ}') + rtMx(nptSel==2,:).*getBetaAll(col{2}, 'acc_{targCrt}');... + (rtMx(nptSel==2,:).^2).*getBetaAll(col{2}, 'acc_{targCrt2}') + (rtMx(nptSel==2,:).^3).*getBetaAll(col{2}, 'acc_{targCrt3}') + (rtMx(nptSel==2,:).^4).*getBetaAll(col{2}, 'acc_{targCrt4}');...
    getBetaAll(col{3}, 'acc_{targ}') + rtMx(nptSel==3,:).*getBetaAll(col{3}, 'acc_{targCrt}');... + (rtMx(nptSel==3,:).^2).*getBetaAll(col{3}, 'acc_{targCrt2}') + (rtMx(nptSel==3,:).^3).*getBetaAll(col{3}, 'acc_{targCrt3}') + (rtMx(nptSel==3,:).^4).*getBetaAll(col{3}, 'acc_{targCrt4}');...
    ];


dyn_dist = (1-dyn_lapse).*[...
    getBetaAll(col{2}, 'acc_{dist}') + rtMx(nptSel==2,:).*getBetaAll(col{2}, 'acc_{distCrt}');... + (rtMx(nptSel==2,:).^2).*getBetaAll(col{2}, 'acc_{distCrt2}') + (rtMx(nptSel==2,:).^3).*getBetaAll(col{2}, 'acc_{distCrt3}') + (rtMx(nptSel==2,:).^4).*getBetaAll(col{2}, 'acc_{distCrt4}');...
    getBetaAll(col{3}, 'acc_{dist}') + rtMx(nptSel==3,:).*getBetaAll(col{3}, 'acc_{distCrt}');... + (rtMx(nptSel==3,:).^2).*getBetaAll(col{3}, 'acc_{distCrt2}') + (rtMx(nptSel==3,:).^3).*getBetaAll(col{3}, 'acc_{distCrt3}') + (rtMx(nptSel==3,:).^4).*getBetaAll(col{3}, 'acc_{distCrt4}');...
    ];



targN = size(dyn_targ,1);
dynN = size(dyn_dist,1);




% plot
n=nexttile([1,10]); hold on;

xline(prctile(acc_rt, 50), '-k');
xline(prctile(acc_rt, 25), '--k');
xline(prctile(acc_rt, 75), '--k');


yline(0, '-k', 'lineWidth', 3);


pltT=plot(rtRange , mean((dyn_targ)), '-g', 'LineWidth', 3);
plot(rtRange , mean((dyn_targ)) + (std((dyn_targ))./sqrt(targN)), '-g', 'LineWidth', .5)
plot(rtRange , mean((dyn_targ)) - (std((dyn_targ))./sqrt(targN)), '-g', 'LineWidth', .5)


pltD=plot(rtRange , mean((dyn_dist)), '-c', 'LineWidth', 3);
plot(rtRange , mean((dyn_dist)) + std((dyn_dist))./sqrt(dynN), '-c', 'LineWidth', .5)
plot(rtRange , mean((dyn_dist)) - std((dyn_dist))./sqrt(dynN), '-c', 'LineWidth', .5)



xlim([min(rtRange), max(rtRange)])
ylim([-1 5]);

xlabel('time (s)')
ylabel('coherence sensitivity')
title('sensitivity dynamics (Attend-Color)')
set(gca, 'TickDir', 'out', 'LineWidth', 1)


















% == plot dynamics LINEAR     ~~~~~~~~~ MOTION ~~~~~~~~~


mdl = mot;
% mdl = exp23;


[rt, crt, d5rt, rt_pred, rt_targ, rt_dist, rt_acc, rt_pt,...
    acc, cacc, acc_pred, acc_targ, acc_dist, acc_rt, acc_pt,...
    all_ndt,...
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
        crt          = [crt; prt - nanmean(prt)];
        
        
        rt_pred     = [rt_pred; out.rt.rt_yhat];
        
        rt_targ     = [acc_targ; discretize(tiedrank(out.data.targ0(out.rt.rt_sel)),nCohDisc)];
        rt_dist     = [acc_dist; discretize(tiedrank(out.data.dist0(out.rt.rt_sel)),nCohDisc)];
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
        
        acc_rt          = [acc_rt; out.data.rt0(out.acc.acc_sel)];
        acc_pt      = [acc_pt; out.acc.acc_X(:,1).*pt_count];
        
        pt_count = pt_count+1;
        
    end
    
    
end



rt_tbl = array2table([rt, crt, rt_pred, rt_targ, rt_dist, rt_acc, rt_pt ], 'VariableNames', ...
    {'rt', 'crt', 'rtY', 'targ', 'dist', 'acc', 'pt'});

acc_tbl = array2table([acc, cacc, acc_pred, acc_rt, d5rt, acc_targ, acc_dist, acc_pt], 'VariableNames', ...
    {'acc', 'cacc', 'accY', 'rt', 'd5rt', 'targ', 'dist', 'pt'});














grpTbl = grpstats(acc_tbl, {'pt'});
% rtMx = repmat([prctile(acc_tbl.rt,5):.001:prctile(acc_tbl.rt,95)], [size(grpTbl,1),1]) - grpTbl.mean_rt;
%
% mRT   = mean(acc_tbl.rt);
% rtRange = (prctile(acc_tbl.rt,5)-mRT):.001:(prctile(acc_tbl.rt,95)-mRT);






rtRange = .4:.001:.6;

clear mRT
for pp = 1:size(grpTbl,1)
    
    mRT(pp,1) = mean(acc_tbl.rt(acc_pt==pp));
    
end
rtMx = repmat(rtRange, [size(grpTbl,1),1]) - mRT;






% estimate
% dyn_targ  = [...
%     getBetaAll(col{1}, 'acc_{targ}') + rtRange.*getBetaAll(col{1}, 'acc_{targCrt}');...
%     getBetaAll(col{2}, 'acc_{targ}') + rtRange.*getBetaAll(col{2}, 'acc_{targCrt}');...
%     getBetaAll(col{3}, 'acc_{targ}') + rtRange.*getBetaAll(col{3}, 'acc_{targCrt}')];
%
% dyn_dist = [...
%     getBetaAll(col{1}, 'acc_{dist}') + rtRange.*getBetaAll(col{1}, 'acc_{distCrt}');...
%     getBetaAll(col{2}, 'acc_{dist}') + rtRange.*getBetaAll(col{2}, 'acc_{distCrt}');...
%     getBetaAll(col{3}, 'acc_{dist}') + rtRange.*getBetaAll(col{3}, 'acc_{distCrt}')];



dyn_lapse  = [...
    1./(1+exp(-(getBetaAll(mot{2}, 'acc_{lapse}') + rtMx(nptSel==2,:).*getBetaAll(mot{2}, 'acc_{crt}'))));
    1./(1+exp(-(getBetaAll(mot{3}, 'acc_{lapse}') + rtMx(nptSel==3,:).*getBetaAll(mot{3}, 'acc_{crt}'))));
    ];




dyn_targ  = (1-dyn_lapse).*[...
    getBetaAll(mot{2}, 'acc_{targ}') + rtMx(nptSel==2,:).*getBetaAll(mot{2}, 'acc_{targCrt}');...
    getBetaAll(mot{3}, 'acc_{targ}') + rtMx(nptSel==3,:).*getBetaAll(mot{3}, 'acc_{targCrt}');...
    ];


dyn_dist = (1-dyn_lapse).*[...
    getBetaAll(mot{2}, 'acc_{dist}') + rtMx(nptSel==2,:).*getBetaAll(mot{2}, 'acc_{distCrt}');...
    getBetaAll(mot{3}, 'acc_{dist}') + rtMx(nptSel==3,:).*getBetaAll(mot{3}, 'acc_{distCrt}');...
    ];



targN = size(dyn_targ,1);
dynN = size(dyn_dist,1);




% plot
n=nexttile([1,10]); hold on;

xline(prctile(acc_rt, 50), '-k');
xline(prctile(acc_rt, 25), '--k');
xline(prctile(acc_rt, 75), '--k');


yline(0, '-k', 'lineWidth', 3);


pltT=plot(rtRange , mean((dyn_targ)), '-g', 'LineWidth', 3);
plot(rtRange , mean((dyn_targ)) + (std((dyn_targ))./sqrt(targN)), '-g', 'LineWidth', .5)
plot(rtRange , mean((dyn_targ)) - (std((dyn_targ))./sqrt(targN)), '-g', 'LineWidth', .5)


pltD=plot(rtRange , mean((dyn_dist)), '-c', 'LineWidth', 3);
plot(rtRange , mean((dyn_dist)) + std((dyn_dist))./sqrt(dynN), '-c', 'LineWidth', .5)
plot(rtRange , mean((dyn_dist)) - std((dyn_dist))./sqrt(dynN), '-c', 'LineWidth', .5)



xlim([min(rtRange), max(rtRange)])
ylim([-2 10]);

xlabel('time (s)')
ylabel('coherence sensitivity')
title('sensitivity dynamics (Attend-Motion)')
set(gca, 'TickDir', 'out', 'LineWidth', 1)
































% % == plot dynamics POLY ==============
% grpTbl = grpstats(acc_tbl, {'pt'});
% % rtMx = repmat([prctile(acc_tbl.rt,5):.001:prctile(acc_tbl.rt,95)], [size(grpTbl,1),1]) - grpTbl.mean_rt;
% %
% % mRT   = mean(acc_tbl.rt);
% % rtRange = (prctile(acc_tbl.rt,5)-mRT):.001:(prctile(acc_tbl.rt,95)-mRT);
%
%
% rtRange = .5:.001:1;
%
% clear mRT
% for pp = 1:size(grpTbl,1)
%
%     mRT(pp,1) = mean(acc_tbl.rt(acc_pt==pp));
%
% end
% rtMx = repmat([.5:.001:1], [size(grpTbl,1),1]) - mRT;
%
%
%
% % estimate
% % dyn_targ  = [...
% %     getBetaAll(col{1}, 'acc_{targ}') + rtRange.*getBetaAll(col{1}, 'acc_{targCrt}');...
% %     getBetaAll(col{2}, 'acc_{targ}') + rtRange.*getBetaAll(col{2}, 'acc_{targCrt}');...
% %     getBetaAll(col{3}, 'acc_{targ}') + rtRange.*getBetaAll(col{3}, 'acc_{targCrt}')];
% %
% % dyn_dist = [...
% %     getBetaAll(col{1}, 'acc_{dist}') + rtRange.*getBetaAll(col{1}, 'acc_{distCrt}');...
% %     getBetaAll(col{2}, 'acc_{dist}') + rtRange.*getBetaAll(col{2}, 'acc_{distCrt}');...
% %     getBetaAll(col{3}, 'acc_{dist}') + rtRange.*getBetaAll(col{3}, 'acc_{distCrt}')];
%
%
%
%
%
%
% dyn_lapse  = [...
%     1./(1+exp(-(getBetaAll(polyMdl{2}, 'acc_{lapse}') + rtMx(nptSel==2,:).*getBetaAll(polyMdl{2}, 'acc_{crt}'))));
%     1./(1+exp(-(getBetaAll(polyMdl{3}, 'acc_{lapse}') + rtMx(nptSel==3,:).*getBetaAll(polyMdl{3}, 'acc_{crt}'))));
%     ];
%
%
%
% dyn_targ  = (1-dyn_lapse).*[...
% %     getBetaAll(polyMdl{1}, 'acc_{targ}') + rtMx(nptSel==1,:).*getBetaAll(polyMdl{1}, 'acc_{targCrt}') + (rtMx(nptSel==1,:).^2).*getBetaAll(polyMdl{1}, 'acc_{targCrt2}') + (rtMx(nptSel==1,:).^3).*getBetaAll(polyMdl{1}, 'acc_{targCrt3}') + (rtMx(nptSel==1,:).^4).*getBetaAll(polyMdl{1}, 'acc_{targCrt4}');...
%     getBetaAll(polyMdl{2}, 'acc_{targ}') + rtMx(nptSel==2,:).*getBetaAll(polyMdl{2}, 'acc_{targCrt}') + (rtMx(nptSel==2,:).^2).*getBetaAll(polyMdl{2}, 'acc_{targCrt2}') + (rtMx(nptSel==2,:).^3).*getBetaAll(polyMdl{2}, 'acc_{targCrt3}') + (rtMx(nptSel==2,:).^4).*getBetaAll(polyMdl{2}, 'acc_{targCrt4}');...
%     getBetaAll(polyMdl{3}, 'acc_{targ}') + rtMx(nptSel==3,:).*getBetaAll(polyMdl{3}, 'acc_{targCrt}') + (rtMx(nptSel==3,:).^2).*getBetaAll(polyMdl{3}, 'acc_{targCrt2}') + (rtMx(nptSel==3,:).^3).*getBetaAll(polyMdl{3}, 'acc_{targCrt3}') + (rtMx(nptSel==3,:).^4).*getBetaAll(polyMdl{3}, 'acc_{targCrt4}');...
%     ];
%
%
% dyn_dist = (1-dyn_lapse).*[...
% %     getBetaAll(polyMdl{1}, 'acc_{dist}') + rtMx(nptSel==1,:).*getBetaAll(polyMdl{1}, 'acc_{distCrt}') + (rtMx(nptSel==1,:).^2).*getBetaAll(polyMdl{1}, 'acc_{distCrt2}') + (rtMx(nptSel==1,:).^3).*getBetaAll(polyMdl{1}, 'acc_{distCrt3}') + (rtMx(nptSel==1,:).^4).*getBetaAll(polyMdl{1}, 'acc_{distCrt4}');...
%     getBetaAll(polyMdl{2}, 'acc_{dist}') + rtMx(nptSel==2,:).*getBetaAll(polyMdl{2}, 'acc_{distCrt}') + (rtMx(nptSel==2,:).^2).*getBetaAll(polyMdl{2}, 'acc_{distCrt2}') + (rtMx(nptSel==2,:).^3).*getBetaAll(polyMdl{2}, 'acc_{distCrt3}') + (rtMx(nptSel==2,:).^4).*getBetaAll(polyMdl{2}, 'acc_{distCrt4}');...
%     getBetaAll(polyMdl{3}, 'acc_{dist}') + rtMx(nptSel==3,:).*getBetaAll(polyMdl{3}, 'acc_{distCrt}') + (rtMx(nptSel==3,:).^2).*getBetaAll(polyMdl{3}, 'acc_{distCrt2}') + (rtMx(nptSel==3,:).^3).*getBetaAll(polyMdl{3}, 'acc_{distCrt3}') + (rtMx(nptSel==3,:).^4).*getBetaAll(polyMdl{3}, 'acc_{distCrt4}');...
%     ];
%
%
%
% targN = size(dyn_targ,1);
% dynN = size(dyn_dist,1);
%
%
%
%
% %
% %
% %
% % pltT=plot(rtRange, mean((dyn_targ)), '--g', 'LineWidth', 3);
% % plot(rtRange, mean((dyn_targ)) + (std((dyn_targ))./sqrt(targN)), '-g', 'LineWidth', .5)
% % plot(rtRange, mean((dyn_targ)) - (std((dyn_targ))./sqrt(targN)), '-g', 'LineWidth', .5)
% %
% %
% % pltD=plot(rtRange, mean((dyn_dist)), '--c', 'LineWidth', 3);
% % plot(rtRange, mean((dyn_dist)) + std((dyn_dist))./sqrt(dynN), '-c', 'LineWidth', .5)
% % plot(rtRange, mean((dyn_dist)) - std((dyn_dist))./sqrt(dynN), '-c', 'LineWidth', .5)
% %
% %
% % xlabel('time (s)')
% % ylabel('coherence sensitivity')
% % title('sensitivity dynamics')
% % set(gca, 'TickDir', 'out', 'LineWidth', 1)
% %
%










% % == plot dynamics EXP ==============
%
%
% mRT   = mean(acc_tbl.rt);
% rtRange = (prctile(acc_tbl.rt,5)):.001:(prctile(acc_tbl.rt,95));
%
%
% dyn_targ  = [...
% %     getBetaAll(expMdl{1}, 'acc_{targFinal}') + (getBetaAll(expMdl{1}, 'acc_{targInit}')-getBetaAll(expMdl{1}, 'acc_{targFinal}')).*exp(-getBetaAll(expMdl{1}, 'acc_{targGain}').*(rtRange));...
%     getBetaAll(expMdl{2}, 'acc_{targFinal}') + (getBetaAll(expMdl{2}, 'acc_{targInit}')-getBetaAll(expMdl{2}, 'acc_{targFinal}')).*exp(-getBetaAll(expMdl{2}, 'acc_{targGain}').*(rtRange));...
%     getBetaAll(expMdl{3}, 'acc_{targFinal}') + (getBetaAll(expMdl{3}, 'acc_{targInit}')-getBetaAll(expMdl{3}, 'acc_{targFinal}')).*exp(-getBetaAll(expMdl{3}, 'acc_{targGain}').*(rtRange));...
%     ];
%
%
% dyn_dist = [...
% %     getBetaAll(expMdl{1}, 'acc_{distFinal}') + (getBetaAll(expMdl{1}, 'acc_{distInit}')-getBetaAll(expMdl{1}, 'acc_{distFinal}')).*exp(-getBetaAll(expMdl{1}, 'acc_{distGain}').*(rtRange));...
%     getBetaAll(expMdl{2}, 'acc_{distFinal}') + (getBetaAll(expMdl{2}, 'acc_{distInit}')-getBetaAll(expMdl{2}, 'acc_{distFinal}')).*exp(-getBetaAll(expMdl{2}, 'acc_{distGain}').*(rtRange));...
%     getBetaAll(expMdl{3}, 'acc_{distFinal}') + (getBetaAll(expMdl{3}, 'acc_{distInit}')-getBetaAll(expMdl{3}, 'acc_{distFinal}')).*exp(-getBetaAll(expMdl{3}, 'acc_{distGain}').*(rtRange));...
%     ];
%
%
%
%
% % dyn_targ  = [...
% %     getBetaAll(expMdl{1}, 'acc_{targFinal}') + (getBetaAll(expMdl{1}, 'acc_{targInit}')-getBetaAll(expMdl{1}, 'acc_{targFinal}')).*exp(-getBetaAll(expMdl{1}, 'acc_{targGain}').*(rtRange-getBetaAll(expMdl{1}, 'acc_{ndt}')));...
% %     getBetaAll(expMdl{2}, 'acc_{targFinal}') + (0-getBetaAll(expMdl{2}, 'acc_{targFinal}')).*exp(-getBetaAll(expMdl{2}, 'acc_{targGain}').*(rtRange-getBetaAll(expMdl{2}, 'acc_{ndt}')));...
% %     getBetaAll(expMdl{3}, 'acc_{targFinal}') + (0-getBetaAll(expMdl{3}, 'acc_{targFinal}')).*exp(-getBetaAll(expMdl{3}, 'acc_{targGain}').*(rtRange-getBetaAll(expMdl{3}, 'acc_{ndt}')));...
% %     ];
% %
% %
% % dyn_dist = [...
% %     0 + (getBetaAll(expMdl{1}, 'acc_{distInit}')-0).*exp(-getBetaAll(expMdl{1}, 'acc_{distGain}').*(rtRange-getBetaAll(expMdl{1}, 'acc_{ndt}')));...
% %     0 + (getBetaAll(expMdl{2}, 'acc_{distInit}')-0).*exp(-getBetaAll(expMdl{2}, 'acc_{distGain}').*(rtRange-getBetaAll(expMdl{2}, 'acc_{ndt}')));...
% %     0 + (getBetaAll(expMdl{3}, 'acc_{distInit}')-0).*exp(-getBetaAll(expMdl{3}, 'acc_{distGain}').*(rtRange-getBetaAll(expMdl{3}, 'acc_{ndt}')));...
% %     ];
%
%
%
%
%
% targN = size(dyn_targ,1);
% dynN = size(dyn_dist,1);
%
%
%
%
% % plot
%
%
% pltT=plot(rtRange, mean((dyn_targ)), '--g', 'LineWidth', 3);
% plot(rtRange, mean((dyn_targ)) + (std((dyn_targ))./sqrt(targN)), '-g', 'LineWidth', .5)
% plot(rtRange, mean((dyn_targ)) - (std((dyn_targ))./sqrt(targN)), '-g', 'LineWidth', .5)
%
%
% pltD=plot(rtRange, mean((dyn_dist)), '--c', 'LineWidth', 3);
% plot(rtRange, mean((dyn_dist)) + std((dyn_dist))./sqrt(dynN), '-c', 'LineWidth', .5)
% plot(rtRange, mean((dyn_dist)) - std((dyn_dist))./sqrt(dynN), '-c', 'LineWidth', .5)
%






% keyboard;



% == plot dynamics EXP 23 ==============


% mRT   = mean(acc_tbl.rt);
% rtRange = (prctile(acc_tbl.rt,5)):.001:(prctile(acc_tbl.rt,95));
%
%

%
% rtRange = .5:.001:1;
%
% clear mRT
% for pp = 1:size(grpTbl,1)
%
%     mRT(pp,1) = mean(acc_tbl.rt(acc_pt==pp));
%
% end
% rtMx = repmat([.5:.001:1], [size(grpTbl,1),1]) - mRT;
%
%
%
%
% dyn_lapse  = [...
%     1./(1+exp(-(getBetaAll(exp23, 'acc_{lapse}') + rtMx(ismember(nptSel, [2,3]),:).*getBetaAll(exp23, 'acc_{crt}'))));
%     ];
%
%
% dyn_targ  = (1-dyn_lapse).*[...
%     getBetaAll(exp23, 'acc_{targFinal}') + (getBetaAll(exp23, 'acc_{targInit}')-getBetaAll(exp23, 'acc_{targFinal}')).*exp(-getBetaAll(exp23, 'acc_{targGain}').*max(0, rtRange - getBetaAll(exp23, 'acc_{targNDT}')));...
%     ];
%
%
% dyn_dist  = (1-dyn_lapse).*[...
%     getBetaAll(exp23, 'acc_{distFinal}') + (getBetaAll(exp23, 'acc_{distInit}')-getBetaAll(exp23, 'acc_{distFinal}')).*exp(-getBetaAll(exp23, 'acc_{distGain}').*max(0, rtRange - getBetaAll(exp23, 'acc_{distNDT}')));...
%     ];
%
%
%
%
% targN = size(dyn_targ,1);
% dynN = size(dyn_dist,1);
%
%
%
%
% % plot
%
%
% pltT=plot(rtRange, mean((dyn_targ)), '-g', 'LineWidth', 1);
% plot(rtRange, mean((dyn_targ)) + (std((dyn_targ))./sqrt(targN)), '-g', 'LineWidth', .5)
% plot(rtRange, mean((dyn_targ)) - (std((dyn_targ))./sqrt(targN)), '-g', 'LineWidth', .5)
%
%
% pltD=plot(rtRange, mean((dyn_dist)), '-c', 'LineWidth', 1);
% plot(rtRange, mean((dyn_dist)) + std((dyn_dist))./sqrt(dynN), '-c', 'LineWidth', .5)
% plot(rtRange, mean((dyn_dist)) - std((dyn_dist))./sqrt(dynN), '-c', 'LineWidth', .5)











% ==== compare fit
%
% for mm = 2:3
%
%     [alpha,exp_r,xp,pxp,bor] = spm_BMS(-[getAIC(col{mm}), getAIC(expMdl{mm}), getAIC(polyMdl{mm})])
%
% end












% === DDM figures





% n=nexttile([1,2]); hold on;
%
%
% yline(2, '-k', 'lineWidth', 1);
% yline(-2, '-k', 'lineWidth', 1);
% yline(0, '--k', 'lineWidth', 1);
%
%
% ddmsim = cumsum(mean(dyn_targ).*log(1.5) - mean(dyn_dist).*1 + randn(size(mean(dyn_dist)))*5)./(sum(mean(dyn_targ).*log(1.5) - mean(dyn_dist).*1)/2);
%
%
% ts = round(linspace(1,length(ddmsim),10));
% dt = mean(diff(ts));
%
% quiver(rtRange(ts) + mRT, ddmsim(ts), mean(dyn_targ(:,ts))*0, (mean(dyn_targ(:,ts)).*log(1.3))/2.5, 0, '-g', 'LineWidth', 1, 'ShowArrowHead', 'on')
% quiver(rtRange(ts) + mRT, ddmsim(ts), mean(dyn_targ(:,ts))*0, -(mean(dyn_dist(:,ts)).*1)/2.5, 0, '-c', 'LineWidth', 1, 'ShowArrowHead', 'on')
% plot(rtRange + mRT, ddmsim, '-k', 'LineWidth', 2);
%
%
% % plot(mean(dyn_targ).*log(1.3), '-g', 'LineWidth', 1)
% % plot(-mean(dyn_dist).*.75, '-c', 'LineWidth', 1)
%
% ylim([-2.2, inf])
% yticks([-2 2])
% yticklabels({'distractor', 'target'})
% % xticks([])
% xlim([min(rtRange+ mRT)-.05, max(rtRange+ mRT)+.05])
%
% xlabel('time')
% set(gca, 'TickDir', 'out', 'LineWidth', 1)


%
%
% n=nexttile([1,8]); hold on;
%
%
% yline(2, '-k', 'lineWidth', 1);
% yline(-2, '-k', 'lineWidth', 1);
% yline(0, '--k', 'lineWidth', 1);
%
%
% ddmsim = cumsum(mean(dyn_targ).*.5 - mean(dyn_dist).*1)./(sum(mean(dyn_targ).*.5 - mean(dyn_dist).*1)/2);
%
%
% ts = round(linspace(1,length(rtRange),9));
% ts = ts(1:end-1);
% dt = mean(diff(rtRange(ts)));
%
% quiver(rtRange(ts) , ddmsim(ts), mean(dyn_targ(:,ts))*0 +(dt), (mean(dyn_targ(:,ts)))/2.5,  '-g', 'LineWidth', 1.5, 'ShowArrowHead', 'off')
% quiver(rtRange(ts) , ddmsim(ts), mean(dyn_targ(:,ts))*0 +(dt), -(mean(dyn_dist(:,ts)).*1)/2.5,  '-c', 'LineWidth', 1.5, 'ShowArrowHead', 'off')
% plot(rtRange , ddmsim, '-k', 'LineWidth', 2);
%
%
% % plot(mean(dyn_targ).*log(1.3), '-g', 'LineWidth', 1)
% % plot(-mean(dyn_dist).*.75, '-c', 'LineWidth', 1)
%
% ylim([-2.2, inf])
% yticks([-2 2])
% yticklabels({'distractor', 'target'})
% % xticks([])
% xlim([min(rtRange), max(rtRange)+.05])
%
% xlabel('time')
% set(gca, 'TickDir', 'out', 'LineWidth', 1)
%
% title('DDM dynamics')






%% ====== CONFLICT ADAPTATION =========================================================
clear col polyMdl expMdl


% tanh AR softplus
col{1} = load('2021-01-18_02-26_exp1_adaptStatic.mat');
col{2} = load('2021-01-20_05-44_exp2_adaptStatic.mat');
col{3} = load('2021-01-18_22-53_exp3_adaptStatic.mat');




mot{1} = load('2021-02-04_19-48_exp1_adaptColMotStatic.mat');
mot{2} = load('2021-02-04_21-45_exp2_adaptColMotStatic.mat');
mot{3} = load('2021-02-05_02-27_exp3_adaptColMotStatic.mat');



dynMdl{1} = load('2021-04-27_17-25_exp1_adaptDualDyn.mat');
dynMdl{2} = load('2021-02-02_06-22_exp2_adaptDualDyn.mat');
dynMdl{3} = load('2021-02-03_02-03_exp3_adaptDualDyn.mat');





disp('');
disp('loaded');
disp('');




% save stats
clc;

disp_mdls = {'col', 'mot' 'dynMdl'}
mdlNames = {'attend-color', 'attend-motion', 'dynamics'};

for cc = 1:3
    
    for mm = 1:length(disp_mdls)
        
        % transform parameters
        x = eval([disp_mdls{mm},'{cc}.results.stats.alpha']);
        for xx = 1:length(x)
           x(xx) =  eval([disp_mdls{mm},'{cc}.data(1).link{xx}(x(xx))']);
        end
        
        
        % correct parameter names
        n = eval([disp_mdls{mm},'{cc}.results.paramName'])';
        for nn = 1:length(n)
            n(nn) = strrep(n(nn), 'acc', 'choice');
            n(nn) = strrep(n(nn), 'rt_{choice}', 'rt_{acc}');
            n(nn) = strrep(n(nn), 'rt_{choiceTarg}', 'rt_{accTarg}');
            n(nn) = strrep(n(nn), 'rt_{choiceDist}', 'rt_{accDist}');
            n(nn) = strrep(n(nn), 'pwr', 'preGain');
            n(nn) = strrep(n(nn), 'dist1', 'prevdist');
            n(nn) = strrep(n(nn), 'Dist1', 'Prevdist');
            n(nn) = strrep(n(nn), 'targ1', 'prevtarg');
            n(nn) = strrep(n(nn), 'Targ1', 'Prevtarg');
            n(nn) = strrep(n(nn), 'crt', 'RT');
        end
        %         disp(n)
        
        % == make table
        tbl = table;
        
        tbl.paramName           = n;
        tbl.df                  = eval(['repmat(length(',disp_mdls{mm},'{cc}.data) - length(', disp_mdls{mm},'{cc}.results.stats.tval), [ length(', disp_mdls{mm},'{cc}.results.stats.tval), 1])']);
        tbl.paramVal            = x;
        tbl.tVal                = eval([disp_mdls{mm},'{cc}.results.stats.tval']);
        tbl.pVal                = eval([disp_mdls{mm},'{cc}.results.stats.p']);
        tbl.cohenD              = eval([disp_mdls{mm},'{cc}.results.stats.alpha./sqrt(diag(', disp_mdls{mm}, '{cc}.results.stats.groupvar))']);
        
        
        
        
        
        
        fprintf('\n =============== exp %d %s =============== \n\n', cc, disp_mdls{mm});
        %         disp(tbl)
        
        writetable(tbl, './tables/adaptation.xlsx', 'Sheet', sprintf('exp%d %s', cc, mdlNames{mm}))
        
    end
    
    %      fprintf('\n\n\n\n\n');
    
end



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

    combP1 = 2*normcdf(nansum(norminv(1-((realmin+pvals)./2)).*df, 2) ./ sqrt(isfinite(pvals)*(df.^2)'), 'upper');
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
clear npt
for ii = 1:3
            npt(cc) = length(col{cc}.data);
end

tcrt = ismember(dynMdl{2}.results.paramName, 'acc_{dist1CrtTarg}');
cv2 = dynMdl{2}.results.stats.groupmeancovariance(tcrt,tcrt);
cv3 = dynMdl{3}.results.stats.groupmeancovariance(tcrt,tcrt);

targ_tval = (dynMdl{2}.results.alpha(tcrt) - dynMdl{3}.results.alpha(tcrt))./...
    sqrt(cv2 + cv3);

targ_df = (cv2 + cv3).^2 ./...
    ((cv2.^2)./(npt(2)-1) + (cv3.^2)./(npt(3)-1));

targ_pval = 2*tcdf(abs(targ_tval), targ_df, 'upper');


%
tcrt = ismember(dynMdl{2}.results.paramName, 'acc_{dist1CrtDist}');
cv2 = dynMdl{2}.results.stats.groupmeancovariance(tcrt,tcrt);
cv3 = dynMdl{3}.results.stats.groupmeancovariance(tcrt,tcrt);

dist_tval = (dynMdl{2}.results.alpha(tcrt) - dynMdl{3}.results.alpha(tcrt))./...
    sqrt(cv2 + cv3);

dist_df = (cv2 + cv3).^2 ./...
    ((cv2.^2)./(npt(2)-1) + (cv3.^2)./(npt(3)-1));

dist_pval = 2*tcdf(abs(dist_tval), dist_df, 'upper');


fprintf('\nexp2&3 targ adapt dyn: t(%.5g)=%.5g, p=%.5g',  targ_df, targ_tval, targ_pval)
fprintf('\nexp2&3 dist adapt dyn: t(%.5g)=%.5g, p=%.5g\n\n',  dist_df, dist_tval, dist_pval)



lg = @(x) 1./(1+exp(-x));

for ii = 1:3
    
    lapse_conInc= 100*[...
        lg(col{ii}.results.alpha(ismember(col{ii}.data(1).paramName, 'acc_{lapse}')) + col{ii}.results.alpha(ismember(col{ii}.data(1).paramName, 'acc_{dist1}'))); ...
        lg(col{ii}.results.alpha(ismember(col{ii}.data(1).paramName, 'acc_{lapse}')) - col{ii}.results.alpha(ismember(col{ii}.data(1).paramName, 'acc_{dist1}'))) ...
       ]
   
end


%% plot conflict adaptation






% load dynamics effects

load('batlow.mat');
bat_cols = batlow;
bat_cIdx = round(linspace(1, length(bat_cols), 7));


mdl = dynMdl;


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

npt2 = [size(col{2}.data,2)];
npt = [size(col{3}.data,2)];

dyn_rt_tbl = array2table([rt, rt_acc, crt, rt_pred, rt_targ, rt_dist, rt_rew, rt_pt ], 'VariableNames', ...
    {'rt', 'acc', 'crt', 'rtY', 'targ', 'dist', 'rew', 'pt'});

dyn_acc_tbl = array2table([acc, cacc, acc_pred, acc_targ, acc_dist, acc_rt, acc_crt, d5rt, acc_rew, acc_dist1, acc_pt], 'VariableNames', ...
    {'acc', 'cacc', 'accY', 'targ', 'dist', 'rt', 'crt', 'drt', 'rew', 'dist1', 'pt'});













mdl = col;


[rt, crt, rt_pred, rt_targ, rt_targ1, rt_dist, rt_dist1, rt_pt,...
    acc_sDist0, acc_sDist1,...
    acc, cacc, acc_pred, acc_targ, acc_targ1, acc_dist, acc_dist1, acc_pt,...
    all_ndt,...
    rt_blk, acc_blk, acc_sBlk, acc_rt, ...
    ] = deal([]);

pt_count = 1;

nCohDisc = 5;




% load('roma.mat');
% cols = roma;
% % cols = colormap('jet');
% cIdx = (round(linspace(1, length(cols), nCohDisc+3)));

cols = [
    231, 111, 81
    244, 162, 97
    233, 196, 106
    42, 157, 143
    38, 70, 83
    ]/255;

cIdx = [0:5];

rtacc_cols = {[36, 123, 160]/255, [242, 95, 92]/255};


whichExps = [1,2,3];

for mm = whichExps
    
    data = mdl{mm}.data;
    results = mdl{mm}.results;
    
    param = results.E';
    ndt = ...
        results.link{ismember(results.paramName, 'rt_{ndt}')} ...
        (param(:, ismember(results.paramName, 'rt_{ndt}')));
    
    all_ndt = [all_ndt; ndt];
    
    for tt = 1:size(data,2)
        
        % simulation parameters
        [~, out]=eval([results.likfun, '(transpose(param(tt,:)), data(tt), [], [], 0, [])']);
        
        
        prt         = log((out.rt.rt(out.rt.rt_sel) - ndt(tt)));
        prt(out.rt.rt(out.rt.rt_sel) <= ndt(tt)) = nan;
        rt          = [rt; prt];
        crt         = [crt; prt - nanmean(prt)];
        
        
        rt_pred     = [rt_pred; out.rt.rt_yhat];
        
        rt_targ     = [rt_targ; discretize(tiedrank(out.data.targ0(out.rt.rt_sel)),nCohDisc)];
        rt_dist     = [rt_dist; discretize(tiedrank(out.data.dist0(out.rt.rt_sel)),nCohDisc)];
        rt_targ1     = [rt_targ1; discretize(tiedrank(out.data.targ1(out.rt.rt_sel)),nCohDisc)];
        rt_dist1     = [rt_dist1; discretize(tiedrank(out.data.dist1(out.rt.rt_sel)),nCohDisc)];
        
        
        %         rt_targ     = [rt_targ; ((out.data.d5Targ0(out.rt.rt_sel)))];
        %         rt_dist     = [rt_dist; ((out.data.d5Dist0(out.rt.rt_sel)))];
        %         rt_targ1     = [rt_targ1; ((out.data.d5Targ1(out.rt.rt_sel)))];
        %         rt_dist1     = [rt_dist1; ((out.data.d5Dist1(out.rt.rt_sel)))];
        
        
        
        rt_pt       = [rt_pt; out.rt.rt_X(:,1).*pt_count];
        
        acc         = [acc; out.data.acc0(out.acc.acc_sel)];
        cacc        = [cacc; center(out.data.acc0(out.acc.acc_sel))];
        
        yhat        = [out.acc.acc_yhat];
        yhat(out.data.corrResp(out.acc.acc_sel)==0) = 1-yhat(out.data.corrResp(out.acc.acc_sel)==0);
        acc_pred    = [acc_pred; yhat];
        
        acc_rt          = [acc_rt; out.data.rt0(out.acc.acc_sel)];
        
        acc_targ    = [acc_targ; discretize(tiedrank(out.data.targ0(out.acc.acc_sel)), nCohDisc)];
        acc_dist    = [acc_dist; discretize(tiedrank((out.data.dist0(out.acc.acc_sel))), nCohDisc)];
        acc_targ1    = [acc_targ1; discretize(tiedrank(out.data.targ1(out.acc.acc_sel)), nCohDisc)];
        acc_dist1    = [acc_dist1; discretize(tiedrank((out.data.dist1(out.acc.acc_sel))), nCohDisc)];
        
        
        %         acc_targ     = [acc_targ; ((out.data.d5Targ0(out.acc.acc_sel)))];
        %         acc_dist     = [acc_dist; ((out.data.d5Dist0(out.acc.acc_sel)))];
        %         acc_targ1     = [acc_targ1; ((out.data.d5Targ1(out.acc.acc_sel)))];
        %         acc_dist1     = [acc_dist1; ((out.data.d5Dist1(out.acc.acc_sel)))];
        
        
        acc_sDist0    = [acc_sDist0; (out.data.dist0(out.acc.acc_sel))];
        acc_sDist1    = [acc_sDist1; (out.data.dist1(out.acc.acc_sel))];
        
        
        acc_pt      = [acc_pt; out.acc.acc_X(:,1).*pt_count];
        
        pt_count = pt_count+1;
        
    end
    
    
end


npt = [size(col{1}.data,2), size(col{2}.data,2), size(col{3}.data,2)];

rt_tbl = array2table([rt, crt, rt_pred, rt_targ, rt_targ1, rt_dist, rt_dist1, rt_pt], 'VariableNames', ...
    {'rt', 'crt', 'rtY', 'targ', 'targ1', 'dist', 'dist1', 'pt'});

acc_tbl = array2table([acc, cacc, acc_pred, acc_targ, acc_targ1, acc_dist, acc_dist1, acc_sDist0, acc_sDist1, acc_rt, acc_pt], 'VariableNames', ...
    {'acc', 'cacc', 'accY', 'targ', 'targ1', 'dist', 'dist1', 'sdist0', 'sdist1', 'rt', 'pt'});



% plot
figure('units','inch','position',[0,0,10,8]); hold on;
tiledlayout(2,5, 'TileSpacing', 'none', 'Padding', 'none');





% =================== Dist1 -> Dist


dotOff = linspace(-.15,.15,nCohDisc);

% ============== acc ==============
grpTbl = grpstats(acc_tbl, {'pt','dist', 'dist1'});
grpTbl = grpstats(grpTbl, {'dist', 'dist1'}, {'mean','sem'});

nexttile([1,2]); hold on;
allY = [];
for ii = 1:nCohDisc
    
    adaptSel = grpTbl.dist1 == ii;
    
    e=errorbar([1:nCohDisc] + dotOff(ii), grpTbl.mean_mean_acc(adaptSel), grpTbl.sem_mean_cacc(adaptSel),...
        'o', 'Color', cols(cIdx(ii+1),:), 'MarkerFaceColor', cols(cIdx(ii+1),:), 'LineWidth', 1);
    
    plot([1:nCohDisc] + dotOff(ii), grpTbl.mean_mean_accY(adaptSel), '-', 'LineWidth', 2, 'Color', cols(cIdx(ii+1),:));
    
    allY = [allY; e.YData(:)];
    
end

title('distractor \rightarrow distractor (Accuracy)')
set(gca, 'TickDir', 'out', 'LineWidth', 1)
yticks(ylim)

xlim([.5, nCohDisc+.5])
xticks([])
% ylabel('Accuracy')
xlabel('Current Distractor Congruence')


% ============== rt ==============
grpTbl = grpstats(rt_tbl, {'pt','dist', 'dist1'}, {'mean','median'});
grpTbl = grpstats(grpTbl, {'dist', 'dist1'}, {'mean','sem'});

nexttile([1,2]); hold on;
allY = [];

for ii = 1:nCohDisc
    
    adaptSel = grpTbl.dist1 == ii;
    
    e=errorbar([1:nCohDisc] + dotOff(ii), exp(grpTbl.mean_mean_rt(adaptSel)) + mean(all_ndt), exp(grpTbl.mean_mean_rt(adaptSel)).*(grpTbl.sem_mean_crt(adaptSel)),...
        'o', 'Color', cols(cIdx(ii+1),:), 'MarkerFaceColor', cols(cIdx(ii+1),:), 'LineWidth', 1);
    
    plot([1:nCohDisc] + dotOff(ii), exp(grpTbl.mean_mean_rtY(adaptSel)) + mean(all_ndt), '-', 'LineWidth', 2, 'Color', cols(cIdx(ii+1),:));
    
    allY = [allY; e.YData(:)];
    
end

set(gca, 'TickDir', 'out', 'LineWidth', 1)
title('distractor \rightarrow distractor (RT)')

yticks(ylim)

xlim([.5, nCohDisc+.5])
xticks([])
xlabel('Current Distractor Congruence')






% ==== betas

x_col_dist_rt  = [1, 1.75, 2.5]+.125;
x_col_dist_acc = [1, 1.75, 2.5]-.125;



M_col_dist_rt = [...
    getBetaMean(col{1}, 'rt_{distDist1}'),    ...
    ...
    getBetaMean(col{2}, 'rt_{distDist1}'),    ...
    ...
    getBetaMean(col{3}, 'rt_{distDist1}'),    ...
    ];

CI_col_dist_rt = [...
    getBetaCI(col{1}, 'rt_{distDist1}'),    ...
    ...
    getBetaCI(col{2}, 'rt_{distDist1}'),    ...
    ...
    getBetaCI(col{3}, 'rt_{distDist1}'),    ...
    ];


all_col_dist_rt{1} = getBetaAll(col{1}, 'rt_{distDist1}');
all_col_dist_rt{2} = getBetaAll(col{2}, 'rt_{distDist1}');
all_col_dist_rt{3} = getBetaAll(col{3}, 'rt_{distDist1}');



M_col_dist_acc = [...
    getBetaMean(col{1}, 'acc_{distDist1}'),    ...
    ...
    getBetaMean(col{2}, 'acc_{distDist1}'),    ...
    ...
    getBetaMean(col{3}, 'acc_{distDist1}'),    ...
    ];

CI_col_dist_acc = [...
    getBetaCI(col{1}, 'acc_{distDist1}'),    ...
    ...
    getBetaCI(col{2}, 'acc_{distDist1}'),    ...
    ...
    getBetaCI(col{3}, 'acc_{distDist1}'),    ...
    ];


all_col_dist_acc{1} = getBetaAll(col{1}, 'acc_{distDist1}');
all_col_dist_acc{2} = getBetaAll(col{2}, 'acc_{distDist1}');
all_col_dist_acc{3} = getBetaAll(col{3}, 'acc_{distDist1}');


n=nexttile([1,1]); hold on;



% acc
yyaxis('left');

for ee = whichExps
    plot(repmat(x_col_dist_acc(ee), [npt(ee),1])+randn([npt(ee),1]).*.01.*ksdensity(all_col_dist_acc{ee},all_col_dist_acc{ee}), all_col_dist_acc{ee},   '.', 'color', rtacc_cols{1});
end

errorbar(x_col_dist_acc(whichExps), M_col_dist_acc(whichExps), CI_col_dist_acc(whichExps),  'ok', 'MarkerFaceColor', 'k', 'LineWidth', 2, 'MarkerSize', 7);

ylim([-.5, 1.5]);
yticks([-.5, 1.5])


% rt
yyaxis('right');
for ee = whichExps
    plot(repmat(x_col_dist_rt(ee), [npt(ee),1])+randn([npt(ee),1]).*.005.*ksdensity(all_col_dist_rt{ee},all_col_dist_rt{ee}), all_col_dist_rt{ee},   '.', 'color', rtacc_cols{2});
end

errorbar(x_col_dist_rt(whichExps), M_col_dist_rt(whichExps), CI_col_dist_rt(whichExps),  'ok', 'MarkerFaceColor', 'k', 'LineWidth', 2, 'MarkerSize', 7);

ylim([-.12, 3*.12]);
yticks([-.12, 3*.12])


% n.YLim = [0, Inf];
set(gca, 'TickDir', 'out', 'LineWidth', 1)
xlim([1-.125-.25, 2.5+.125+.25])
yline(0, '-k', 'LineWidth', 1);

xticks([1, 1.75, 2.5])
xticklabels({'Exp1', 'Exp2', 'Exp3'})










% ============== Dist -> Targ ==============

% ============== D -> T acc ==============
grpTbl = grpstats(acc_tbl, {'pt','targ', 'dist1'});
grpTbl = grpstats(grpTbl, {'targ', 'dist1'}, {'mean','sem'});

nexttile([1,2]); hold on;

allY = [];
for ii = 1:nCohDisc
    
    adaptSel = grpTbl.dist1 == ii;
    
    e=errorbar([1:nCohDisc] + dotOff(ii), grpTbl.mean_mean_acc(adaptSel), grpTbl.sem_mean_cacc(adaptSel),...
        'o', 'Color', cols(cIdx(ii+1),:), 'MarkerFaceColor', cols(cIdx(ii+1),:), 'LineWidth', 1);
    
    plot([1:nCohDisc] + dotOff(ii), grpTbl.mean_mean_accY(adaptSel), '-', 'LineWidth', 2, 'Color', cols(cIdx(ii+1),:));
    
    allY = [allY; e.YData(:)];
end

yticks(ylim)
title('distractor \rightarrow target (Accuracy)')
set(gca, 'TickDir', 'out', 'LineWidth', 1)
xlim([.5, nCohDisc+.5])
xticks([])
xlabel('Current Target Coherence')









% ============== D -> T RT ==============
grpTbl = grpstats(rt_tbl, {'pt','targ', 'dist1'});
grpTbl = grpstats(grpTbl, {'targ', 'dist1'}, {'mean','sem'});

nexttile([1,2]); hold on;

allY = [];
for ii = 1:nCohDisc
    
    adaptSel = grpTbl.dist1 == ii;
    
    e=errorbar([1:nCohDisc] + dotOff(ii), exp(grpTbl.mean_mean_rt(adaptSel)) + mean(all_ndt), exp(grpTbl.mean_mean_rt(adaptSel)).*(grpTbl.sem_mean_crt(adaptSel)),...
        'o', 'Color', cols(cIdx(ii+1),:), 'MarkerFaceColor', cols(cIdx(ii+1),:), 'LineWidth', 1);
    
    plot([1:nCohDisc] + dotOff(ii), exp(grpTbl.mean_mean_rtY(adaptSel)) + mean(all_ndt), '-', 'LineWidth', 2, 'Color', cols(cIdx(ii+1),:));
    
    allY = [allY; e.YData(:)];
end

yticks(ylim)
title('distractor \rightarrow target (RT)')
set(gca, 'TickDir', 'out', 'LineWidth', 1)
xlim([.5, nCohDisc+.5])
xticks([])
xlabel('Current Target Coherence')





% plot betas

x_col_dist_rt  = [1, 1.75, 2.5]+.125;
x_col_dist_acc = [1, 1.75, 2.5]-.125;



M_col_dist_rt = [...
%     getBetaMean(col{1}, 'rt_{targDist1}'),    ...
    ...
    getBetaMean(col{2}, 'rt_{targDist1}'),    ...
    ...
    getBetaMean(col{3}, 'rt_{targDist1}'),    ...
    ];

CI_col_dist_rt = [...
%     getBetaCI(col{1}, 'rt_{targDist1}'),    ...
    ...
    getBetaCI(col{2}, 'rt_{targDist1}'),    ...
    ...
    getBetaCI(col{3}, 'rt_{targDist1}'),    ...
    ];


% all_col_dist_rt{1} = getBetaAll(col{1}, 'rt_{targDist1}');
all_col_dist_rt{2} = getBetaAll(col{2}, 'rt_{targDist1}');
all_col_dist_rt{3} = getBetaAll(col{3}, 'rt_{targDist1}');



M_col_dist_acc = [...
%     getBetaMean(col{1}, 'acc_{targDist1}'),    ...
    ...
    getBetaMean(col{2}, 'acc_{targDist1}'),    ...
    ...
    getBetaMean(col{3}, 'acc_{targDist1}'),    ...
    ];

CI_col_dist_acc = [...
%     getBetaCI(col{1}, 'acc_{targDist1}'),    ...
    ...
    getBetaCI(col{2}, 'acc_{targDist1}'),    ...
    ...
    getBetaCI(col{3}, 'acc_{targDist1}'),    ...
    ];


% all_col_dist_acc{1} = getBetaAll(col{1}, 'acc_{targDist1}');
all_col_dist_acc{2} = getBetaAll(col{2}, 'acc_{targDist1}');
all_col_dist_acc{3} = getBetaAll(col{3}, 'acc_{targDist1}');


n=nexttile([1,1]); hold on;



% acc
yyaxis('left');
whichExps = [2:3];
for ee = whichExps
    plot(repmat(x_col_dist_acc(ee), [npt(ee),1])+randn([npt(ee),1]).*.01.*ksdensity(all_col_dist_acc{ee},all_col_dist_acc{ee}), all_col_dist_acc{ee},   '.', 'color', rtacc_cols{1});
end

errorbar(x_col_dist_acc(whichExps), M_col_dist_acc, CI_col_dist_acc,  'ok', 'MarkerFaceColor', 'k', 'LineWidth', 2, 'MarkerSize', 7);

ylim([-.4, .4]);
yticks(ylim)
xlim([.5, nCohDisc+.5])
xticks([])
xlabel('Current Target Coherence')

% rt
yyaxis('right');
for ee = whichExps
    plot(repmat(x_col_dist_rt(ee), [npt(ee),1])+randn([npt(ee),1]).*.005.*ksdensity(all_col_dist_rt{ee},all_col_dist_rt{ee}), all_col_dist_rt{ee},   '.', 'color', rtacc_cols{2});
end

errorbar(x_col_dist_rt(whichExps), M_col_dist_rt, CI_col_dist_rt,  'ok', 'MarkerFaceColor', 'k', 'LineWidth', 2, 'MarkerSize', 7);

ylim([-.2, .2]);
yticks(ylim)
xlim([.5, nCohDisc+.5])
xticks([])
xlabel('Current Target Coherence')

% n.YLim = [0, Inf];
set(gca, 'TickDir', 'out', 'LineWidth', 1)
xlim([1-.125-.25, 2.5+.125+.25])
yline(0, '-k', 'LineWidth', 1);

xticks([1, 1.75, 2.5])
xticklabels({'Exp1', 'Exp2', 'Exp3'})










% plot 2

figure('units','inch','position',[0,0,10,8]); hold on;
tiledlayout(2,5, 'TileSpacing', 'none', 'Padding', 'none');


% ============== T -> T acc ==============
grpTbl = grpstats(acc_tbl, {'pt','targ', 'targ1'});
grpTbl = grpstats(grpTbl, {'targ', 'targ1'}, {'mean','sem'});

nexttile([1,2]); hold on;
allY = [];
for ii = 1:nCohDisc
    
    adaptSel = grpTbl.targ1 == ii;
    
    e=errorbar([1:nCohDisc] + dotOff(ii), grpTbl.mean_mean_acc(adaptSel), grpTbl.sem_mean_cacc(adaptSel),...
        'o', 'Color', cols(cIdx(ii+1),:), 'MarkerFaceColor', cols(cIdx(ii+1),:), 'LineWidth', 1);
    
    plot([1:nCohDisc] + dotOff(ii), grpTbl.mean_mean_accY(adaptSel), '-', 'LineWidth', 2, 'Color', cols(cIdx(ii+1),:));
    
    allY = [allY; e.YData(:)];
    
end

title('target \rightarrow target (Accuracy)')
set(gca, 'TickDir', 'out', 'LineWidth', 1)
yticks(ylim)
xlim([.5, nCohDisc+.5])
xticks([])
xlabel('Current Target Coherence')









% ============== T -> T RT ==============
grpTbl = grpstats(rt_tbl, {'pt','targ', 'targ1'});
grpTbl = grpstats(grpTbl, {'targ', 'targ1'}, {'mean','sem'});

nexttile([1,2]); hold on;
allY = [];
for ii = 1:nCohDisc
    
    adaptSel = grpTbl.targ1 == ii;
    
    e=errorbar([1:nCohDisc] + dotOff(ii), exp(grpTbl.mean_mean_rt(adaptSel)) + mean(all_ndt), exp(grpTbl.mean_mean_rt(adaptSel)).*(grpTbl.sem_mean_crt(adaptSel)),...
        'o', 'Color', cols(cIdx(ii+1),:), 'MarkerFaceColor', cols(cIdx(ii+1),:), 'LineWidth', 1);
    
    plot([1:nCohDisc] + dotOff(ii), exp(grpTbl.mean_mean_rtY(adaptSel)) + mean(all_ndt), '-', 'LineWidth', 2, 'Color', cols(cIdx(ii+1),:));
    
    allY = [allY; e.YData(:)];
    
end

title('target \rightarrow target (RT)')
set(gca, 'TickDir', 'out', 'LineWidth', 1)
yticks(ylim)
xlim([.5, nCohDisc+.5])
xticks([])
xlabel('Current Target Coherence')







% plot betas

x_col_dist_rt  = [1, 1.75, 2.5]+.125;
x_col_dist_acc = [1, 1.75, 2.5]-.125;



M_col_dist_rt = [...
%     getBetaMean(col{1}, 'rt_{targDist1}'),    ...
    ...
    getBetaMean(col{2}, 'rt_{targTarg1}'),    ...
    ...
    getBetaMean(col{3}, 'rt_{targTarg1}'),    ...
    ];

CI_col_dist_rt = [...
%     getBetaCI(col{1}, 'rt_{targDist1}'),    ...
    ...
    getBetaCI(col{2}, 'rt_{targTarg1}'),    ...
    ...
    getBetaCI(col{3}, 'rt_{targTarg1}'),    ...
    ];


% all_col_dist_rt{1} = getBetaAll(col{1}, 'rt_{targDist1}');
all_col_dist_rt{2} = getBetaAll(col{2}, 'rt_{targTarg1}');
all_col_dist_rt{3} = getBetaAll(col{3}, 'rt_{targTarg1}');



M_col_dist_acc = [...
%     getBetaMean(col{1}, 'acc_{targDist1}'),    ...
    ...
    getBetaMean(col{2}, 'acc_{targTarg1}'),    ...
    ...
    getBetaMean(col{3}, 'acc_{targTarg1}'),    ...
    ];

CI_col_dist_acc = [...
%     getBetaCI(col{1}, 'acc_{targDist1}'),    ...
    ...
    getBetaCI(col{2}, 'acc_{targTarg1}'),    ...
    ...
    getBetaCI(col{3}, 'acc_{targTarg1}'),    ...
    ];


% all_col_dist_acc{1} = getBetaAll(col{1}, 'acc_{targDist1}');
all_col_dist_acc{2} = getBetaAll(col{2}, 'acc_{targTarg1}');
all_col_dist_acc{3} = getBetaAll(col{3}, 'acc_{targTarg1}');


n=nexttile([1,1]); hold on;



% acc
yyaxis('left');
whichExps = [2:3];
for ee = whichExps
    plot(repmat(x_col_dist_acc(ee), [npt(ee),1])+randn([npt(ee),1]).*.01.*ksdensity(all_col_dist_acc{ee},all_col_dist_acc{ee}), all_col_dist_acc{ee},   '.', 'color', rtacc_cols{1});
end

errorbar(x_col_dist_acc(whichExps), M_col_dist_acc, CI_col_dist_acc,  'ok', 'MarkerFaceColor', 'k', 'LineWidth', 2, 'MarkerSize', 7);

ylim([-2, 2]);
yticks(ylim)


% rt
yyaxis('right');
for ee = whichExps
    plot(repmat(x_col_dist_rt(ee), [npt(ee),1])+randn([npt(ee),1]).*.005.*ksdensity(all_col_dist_rt{ee},all_col_dist_rt{ee}), all_col_dist_rt{ee},   '.', 'color', rtacc_cols{2});
end

errorbar(x_col_dist_rt(whichExps), M_col_dist_rt, CI_col_dist_rt,  'ok', 'MarkerFaceColor', 'k', 'LineWidth', 2, 'MarkerSize', 7);

ylim([-.8, .8]);
yticks(ylim)


% n.YLim = [0, Inf];
set(gca, 'TickDir', 'out', 'LineWidth', 1)
xlim([1-.125-.25, 2.5+.125+.25])
yline(0, '-k', 'LineWidth', 1);

xticks([1, 1.75, 2.5])
xticklabels({'Exp1', 'Exp2', 'Exp3'})














% ============== T -> D acc ==============
grpTbl = grpstats(acc_tbl, {'pt','dist', 'targ1'});
grpTbl = grpstats(grpTbl, {'dist', 'targ1'}, {'mean','sem'});

nexttile([1,2]); hold on;
allY = [];
for ii = 1:nCohDisc
    
    adaptSel = grpTbl.targ1 == ii;
    
    e=errorbar([1:nCohDisc] + dotOff(ii), grpTbl.mean_mean_acc(adaptSel), grpTbl.sem_mean_cacc(adaptSel),...
        'o', 'Color', cols(cIdx(ii+1),:), 'MarkerFaceColor', cols(cIdx(ii+1),:), 'LineWidth', 1);
    
    plot([1:nCohDisc] + dotOff(ii), grpTbl.mean_mean_accY(adaptSel), '-', 'LineWidth', 2, 'Color', cols(cIdx(ii+1),:));
    
    allY = [allY; e.YData(:)];
    
end

title('target \rightarrow distractor (Accuracy)')
set(gca, 'TickDir', 'out', 'LineWidth', 1)
yticks(ylim)
xlim([.5, nCohDisc+.5])
xticks([])
xlabel('Current Distractor Congruence')









% ============== T -> D RT ==============
grpTbl = grpstats(rt_tbl, {'pt','dist', 'targ1'});
grpTbl = grpstats(grpTbl, {'dist', 'targ1'}, {'mean','sem'});

nexttile([1,2]); hold on;
allY = [];
for ii = 1:nCohDisc
    
    adaptSel = grpTbl.targ1 == ii;
    
    e=errorbar([1:nCohDisc] + dotOff(ii), exp(grpTbl.mean_mean_rt(adaptSel)) + mean(all_ndt), exp(grpTbl.mean_mean_rt(adaptSel)).*(grpTbl.sem_mean_crt(adaptSel)),...
        'o', 'Color', cols(cIdx(ii+1),:), 'MarkerFaceColor', cols(cIdx(ii+1),:), 'LineWidth', 1);
    
    plot([1:nCohDisc] + dotOff(ii), exp(grpTbl.mean_mean_rtY(adaptSel)) + mean(all_ndt), '-', 'LineWidth', 2, 'Color', cols(cIdx(ii+1),:));
    
    allY = [allY; e.YData(:)];
    
end

title('target \rightarrow distractor (RT)')
set(gca, 'TickDir', 'out', 'LineWidth', 1)
yticks(ylim)
xlim([.5, nCohDisc+.5])
xticks([])
xlabel('Current Distractor Congruence')












% plot betas

x_col_dist_rt  = [1, 1.75, 2.5]+.125;
x_col_dist_acc = [1, 1.75, 2.5]-.125;



M_col_dist_rt = [...
%     getBetaMean(col{1}, 'rt_{targDist1}'),    ...
    ...
    getBetaMean(col{2}, 'rt_{distTarg1}'),    ...
    ...
    getBetaMean(col{3}, 'rt_{distTarg1}'),    ...
    ];

CI_col_dist_rt = [...
%     getBetaCI(col{1}, 'rt_{targDist1}'),    ...
    ...
    getBetaCI(col{2}, 'rt_{distTarg1}'),    ...
    ...
    getBetaCI(col{3}, 'rt_{distTarg1}'),    ...
    ];


% all_col_dist_rt{1} = getBetaAll(col{1}, 'rt_{targDist1}');
all_col_dist_rt{2} = getBetaAll(col{2}, 'rt_{distTarg1}');
all_col_dist_rt{3} = getBetaAll(col{3}, 'rt_{distTarg1}');



M_col_dist_acc = [...
%     getBetaMean(col{1}, 'acc_{targDist1}'),    ...
    ...
    getBetaMean(col{2}, 'acc_{distTarg1}'),    ...
    ...
    getBetaMean(col{3}, 'acc_{distTarg1}'),    ...
    ];

CI_col_dist_acc = [...
%     getBetaCI(col{1}, 'acc_{targDist1}'),    ...
    ...
    getBetaCI(col{2}, 'acc_{distTarg1}'),    ...
    ...
    getBetaCI(col{3}, 'acc_{distTarg1}'),    ...
    ];


% all_col_dist_acc{1} = getBetaAll(col{1}, 'acc_{targDist1}');
all_col_dist_acc{2} = getBetaAll(col{2}, 'acc_{distTarg1}');
all_col_dist_acc{3} = getBetaAll(col{3}, 'acc_{distTarg1}');


n=nexttile([1,1]); hold on;



% acc
yyaxis('left');
whichExps = [2:3]
for ee = whichExps
    plot(repmat(x_col_dist_acc(ee), [npt(ee),1])+randn([npt(ee),1]).*.01.*ksdensity(all_col_dist_acc{ee},all_col_dist_acc{ee}), all_col_dist_acc{ee},   '.', 'color', rtacc_cols{1});
end

errorbar(x_col_dist_acc(whichExps), M_col_dist_acc, CI_col_dist_acc,  'ok', 'MarkerFaceColor', 'k', 'LineWidth', 2, 'MarkerSize', 7);

ylim([-2, 2]);
yticks(ylim)


% rt
yyaxis('right');
for ee = whichExps
    plot(repmat(x_col_dist_rt(ee), [npt(ee),1])+randn([npt(ee),1]).*.005.*ksdensity(all_col_dist_rt{ee},all_col_dist_rt{ee}), all_col_dist_rt{ee},   '.', 'color', rtacc_cols{2});
end

errorbar(x_col_dist_rt(whichExps), M_col_dist_rt, CI_col_dist_rt,  'ok', 'MarkerFaceColor', 'k', 'LineWidth', 2, 'MarkerSize', 7);

ylim([-.2, .2]);
yticks(ylim)


% n.YLim = [0, Inf];
set(gca, 'TickDir', 'out', 'LineWidth', 1)
xlim([1-.125-.25, 2.5+.125+.25])
yline(0, '-k', 'LineWidth', 1);

xticks([1, 1.75, 2.5])
xticklabels({'Exp1', 'Exp2', 'Exp3'})


%% ============== plot ADAPT Dynamics effects ==============

figure('units','inch','position',[0,0,1.5*7.5,1.5*2]); hold on;
tiledlayout(1,5, 'TileSpacing', 'none', 'Padding', 'none');



grpTbl = grpstats(dyn_acc_tbl, {'pt', 'dist', 'drt' 'dist1'});
grpTbl = grpstats(grpTbl, {'dist', 'drt', 'dist1'}, {'mean','sem'});


nexttile([1,2]); hold on;
yline(0, '--k', 'LineWidth', 1)

xxoff = linspace(-.1,.1,5);
for ii = [1:5]
    
    plot([1:5]+xxoff(ii), (grpTbl.mean_mean_acc(grpTbl.drt == ii & grpTbl.dist1==3))-(grpTbl.mean_mean_acc(grpTbl.drt == ii & grpTbl.dist1==1)),...
        'o', 'LineWidth', 1, 'color', bat_cols(bat_cIdx(ii+1),:),  'MarkerFaceColor', bat_cols(bat_cIdx(ii+1),:));
    
    plot([1:5]+xxoff(ii), (grpTbl.mean_mean_accY(grpTbl.drt == ii & grpTbl.dist1==3))-(grpTbl.mean_mean_accY(grpTbl.drt == ii & grpTbl.dist1==1)),...
        '-', 'LineWidth', 2, 'color', bat_cols(bat_cIdx(ii+1),:));
    
end

set(gca, 'TickDir', 'out', 'LineWidth', 1);
xlim([.5, 5.5])
title('adaptation dynamics')
ylabel('Cong - Inc (accuracy)')
xlabel('distractor congruence')
% legend(pltp, {'fast','','med','','slow'}, 'Location', 'northeast')















% ==== plot dynamics
nptSel = [];
for ww = whichExps
    nptSel = [nptSel; ones(npt(ww),1)*ww];
end
grpTbl = grpstats(acc_tbl, {'pt'});
rtMx = repmat([prctile(acc_tbl.rt,5):.001:1], [size(grpTbl,1),1]) - grpTbl.mean_rt;


mRT   = mean(acc_tbl.rt);
rtRange = [(prctile(acc_tbl.rt,5)):.001:1]-mRT;


% plot
n=nexttile([1,2]); hold on;


xline(prctile(acc_rt, 50), '-k');
xline(prctile(acc_rt, 25), '--k');
xline(prctile(acc_rt, 75), '--k');



% for bb = 1:nCohDisc
%
%     mBlk = mean(acc_sDist1(acc_dist1==bb));
%
%     dyn_targ  = [...
%         getBetaAll(polyMdl{2}, 'acc_{targ}') + rtRange.*getBetaAll(polyMdl{2}, 'acc_{crtTarg}') + mBlk.*getBetaAll(polyMdl{2}, 'acc_{dist1Targ}') + rtRange.*mBlk.*getBetaAll(polyMdl{2}, 'acc_{dist1CrtTarg}');...
%         getBetaAll(polyMdl{3}, 'acc_{targ}') + rtRange.*getBetaAll(polyMdl{3}, 'acc_{crtTarg}') + mBlk.*getBetaAll(polyMdl{3}, 'acc_{dist1Targ}') + rtRange.*mBlk.*getBetaAll(polyMdl{3}, 'acc_{dist1CrtTarg}');...
%         ];
%
%
%    dyn_dist  = [...
%        getBetaAll(polyMdl{1}, 'acc_{dist}') + rtRange.*getBetaAll(polyMdl{1}, 'acc_{crtDist}') + mBlk.*getBetaAll(polyMdl{1}, 'acc_{dist1Dist}') + rtRange.*mBlk.*getBetaAll(polyMdl{1}, 'acc_{dist1CrtDist}');...
%         getBetaAll(polyMdl{2}, 'acc_{dist}') + rtRange.*getBetaAll(polyMdl{2}, 'acc_{crtDist}') + mBlk.*getBetaAll(polyMdl{2}, 'acc_{dist1Dist}') + rtRange.*mBlk.*getBetaAll(polyMdl{2}, 'acc_{dist1CrtDist}');...
%         getBetaAll(polyMdl{3}, 'acc_{dist}') + rtRange.*getBetaAll(polyMdl{3}, 'acc_{crtDist}') + mBlk.*getBetaAll(polyMdl{3}, 'acc_{dist1Dist}') + rtRange.*mBlk.*getBetaAll(polyMdl{3}, 'acc_{dist1CrtDist}');...
%         ];
%
%
%
%     targN = size(dyn_targ,1);
%     dynN = size(dyn_dist,1);
%
%
%
%
%
%     yline(0, '-k', 'lineWidth', 3);
%
%
%     pltT=plot(rtRange + mRT, mean((dyn_targ)), '-', 'LineWidth', 3, 'Color', cols(cIdx(bb+1),:));
%     plot(rtRange + mRT, mean((dyn_targ)) + (std((dyn_targ))./sqrt(targN)), '-', 'Color', cols(cIdx(bb+1),:), 'LineWidth', .5)
%     plot(rtRange + mRT, mean((dyn_targ)) - (std((dyn_targ))./sqrt(targN)), '-', 'Color', cols(cIdx(bb+1),:), 'LineWidth', .5)
%
%
%     pltD=plot(rtRange + mRT, mean((dyn_dist)), '-', 'LineWidth', 3, 'Color', cols(cIdx(bb+1),:));
%     plot(rtRange + mRT, mean((dyn_dist)) + std((dyn_dist))./sqrt(dynN), '-', 'Color', cols(cIdx(bb+1),:),  'LineWidth', .5)
%     plot(rtRange + mRT, mean((dyn_dist)) - std((dyn_dist))./sqrt(dynN), '-', 'Color', cols(cIdx(bb+1),:), 'LineWidth', .5)
%
%
% end




rtRange = .5:.001:1;

colDark = linspace(1,.2,nCohDisc);


for bb = [1:nCohDisc]
    
    mBlk = mean(acc_sDist1(acc_dist1==bb));
    
    clear mRT
    for pp = 1:size(grpTbl,1)
        
%         mRT(pp,1) = (acc_sDist1(acc_pt==pp)\center(acc_tbl.rt(acc_pt==pp)))*mBlk + mean(acc_tbl.rt(acc_pt==pp));
        mRT(pp,1) = mean(dyn_acc_tbl.crt(acc_pt==pp & acc_dist1==bb)) + mean(acc_tbl.rt(acc_pt==pp));
        
    end
    rtMx = repmat(rtRange, [size(grpTbl,1),1]) - mRT;
    
    
    
    
    
    
    
    dyn_lapse  = [...
        1./(1+exp(-(getBetaAll(dynMdl{2}, 'acc_{lapse}') + rtMx(nptSel==2,:).*getBetaAll(dynMdl{2}, 'acc_{crt}') + mBlk.*getBetaAll(dynMdl{2}, 'acc_{dist1}') + mBlk.*rtMx(nptSel==2,:).*getBetaAll(dynMdl{2}, 'acc_{dist1Crt}'))));
        1./(1+exp(-(getBetaAll(dynMdl{3}, 'acc_{lapse}') + rtMx(nptSel==3,:).*getBetaAll(dynMdl{3}, 'acc_{crt}') + mBlk.*getBetaAll(dynMdl{3}, 'acc_{dist1}') + mBlk.*rtMx(nptSel==3,:).*getBetaAll(dynMdl{3}, 'acc_{dist1Crt}'))));
        ];
    
    
    dyn_targ  = (1-dyn_lapse).*[...
        getBetaAll(dynMdl{2}, 'acc_{targ}') + rtMx(nptSel==2,:).*getBetaAll(dynMdl{2}, 'acc_{crtTarg}') + mBlk.*getBetaAll(dynMdl{2}, 'acc_{dist1Targ}') + mBlk.*rtMx(nptSel==2,:).*getBetaAll(dynMdl{2}, 'acc_{dist1CrtTarg}');...
        getBetaAll(dynMdl{3}, 'acc_{targ}') + rtMx(nptSel==3,:).*getBetaAll(dynMdl{3}, 'acc_{crtTarg}') + mBlk.*getBetaAll(dynMdl{3}, 'acc_{dist1Targ}') + mBlk.*rtMx(nptSel==3,:).*getBetaAll(dynMdl{3}, 'acc_{dist1CrtTarg}');...
        ];
    
    
    dyn_dist = (1-dyn_lapse).*[...
        getBetaAll(dynMdl{2}, 'acc_{dist}') + rtMx(nptSel==2,:).*getBetaAll(dynMdl{2}, 'acc_{crtDist}') + mBlk.*getBetaAll(dynMdl{2}, 'acc_{dist1Dist}') + mBlk.*rtMx(nptSel==2,:).*getBetaAll(dynMdl{2}, 'acc_{dist1CrtDist}');...
        getBetaAll(dynMdl{3}, 'acc_{dist}') + rtMx(nptSel==3,:).*getBetaAll(dynMdl{3}, 'acc_{crtDist}') + mBlk.*getBetaAll(dynMdl{3}, 'acc_{dist1Dist}') + mBlk.*rtMx(nptSel==3,:).*getBetaAll(dynMdl{3}, 'acc_{dist1CrtDist}');...
        ];
    
    
    
    
    
    targN = size(dyn_targ,1);
    dynN = size(dyn_dist,1);
    
    
    
    
    
    yline(0, '--k', 'lineWidth', 3);
    
    
    %     pltT=plot(rtRange, mean((dyn_targ)), '-', 'LineWidth', 3, 'Color', cols(cIdx(bb+1),:));
    %     plot(rtRange, mean((dyn_targ)) + (std((dyn_targ))./sqrt(targN)), '-', 'Color', cols(cIdx(bb+1),:), 'LineWidth', .5)
    %     plot(rtRange, mean((dyn_targ)) - (std((dyn_targ))./sqrt(targN)), '-', 'Color', cols(cIdx(bb+1),:), 'LineWidth', .5)
    
    pltT=plot(rtRange, mean((dyn_targ)), '-', 'LineWidth', 3, 'Color', [0, colDark(bb), 0]);
    plot(rtRange, mean((dyn_targ)) + (std((dyn_targ))./sqrt(targN)), '-', 'Color', [0, colDark(bb), 0], 'LineWidth', .5)
    plot(rtRange, mean((dyn_targ)) - (std((dyn_targ))./sqrt(targN)), '-', 'Color', [0, colDark(bb), 0], 'LineWidth', .5)
    
    
    %
    %     pltD=plot(rtRange, mean((dyn_dist)), '-', 'LineWidth', 3, 'Color', cols(cIdx(bb+1),:));
    %     plot(rtRange, mean((dyn_dist)) + std((dyn_dist))./sqrt(dynN), '-', 'Color', cols(cIdx(bb+1),:),  'LineWidth', .5)
    %     plot(rtRange, mean((dyn_dist)) - std((dyn_dist))./sqrt(dynN), '-', 'Color', cols(cIdx(bb+1),:), 'LineWidth', .5)
    
    
    pltD=plot(rtRange, mean((dyn_dist)), '-', 'LineWidth', 3, 'Color', [0, colDark(bb), colDark(bb)]);
    plot(rtRange, mean((dyn_dist)) + std((dyn_dist))./sqrt(dynN), '-', 'Color', [0, colDark(bb), colDark(bb)],  'LineWidth', .5)
    plot(rtRange, mean((dyn_dist)) - std((dyn_dist))./sqrt(dynN), '-', 'Color', [0, colDark(bb), colDark(bb)], 'LineWidth', .5)
    
    
end




xlim([.5, 1])
ylim([-1 6]);
yticks([-1, 6])

xlabel('Time (s)')
ylabel('Coherence Sensitivity')
title('Sensitivity Dynamics')
set(gca, 'TickDir', 'out', 'LineWidth', 1)

% legend([pltT, pltD], {'target', 'distractor'}, 'Location', 'best')






















% ============== targ dynamics betas ==============
x_col_targ_acc = [1.5, 2.25];
x_col_dist_acc = [1, 1.75, 2.5]-.1;


M_col_dist_acc = [...
    %     getBetaMean(dynMdl{1}, 'acc_{dist1CrtTarg}'),    ...
    ...
    getBetaMean(dynMdl{2}, 'acc_{dist1CrtTarg}'),    ...
    ...
    getBetaMean(dynMdl{3}, 'acc_{dist1CrtTarg}'),    ...
    ];

CI_col_dist_acc = [...
    %     getBetaCI(dynMdl{1}, 'acc_{dist1CrtTarg}'),    ...
    ...
    getBetaCI(dynMdl{2}, 'acc_{dist1CrtTarg}'),    ...
    ...
    getBetaCI(dynMdl{3}, 'acc_{dist1CrtTarg}'),    ...
    ];


% all_col_dist_acc{1} = getBetaAll(dynMdl{1}, 'acc_{dist1CrtTarg}')*nan;
all_col_dist_acc{2} = getBetaAll(dynMdl{2}, 'acc_{dist1CrtTarg}');
all_col_dist_acc{3} = getBetaAll(dynMdl{3}, 'acc_{dist1CrtTarg}');


n=nexttile([1,1]); hold on;

for ee = 2:3
    plot(repmat(x_col_dist_acc(ee), [npt(ee),1])+randn([npt(ee),1]).*.1.*ksdensity(all_col_dist_acc{ee},all_col_dist_acc{ee}), all_col_dist_acc{ee},   '.g');
end

errorbar(x_col_dist_acc(2:3), M_col_dist_acc(1:2), CI_col_dist_acc(1:2),  'ok', 'MarkerFaceColor', 'g', 'LineWidth', 2, 'MarkerSize', 7);


% n.YLim = [0, Inf];
set(gca, 'TickDir', 'out', 'LineWidth', 1)
xlim([x_col_dist_acc(1)-.25, x_col_dist_acc(end) + .25])
yline(0, '-k', 'LineWidth', 1);





% ============== dist dynamics betas ==============
x_col_targ_acc = [1.5, 2.25];
x_col_dist_acc = [1, 1.75, 2.5]+.1;


M_col_dist_acc = [...
    getBetaMean(dynMdl{1}, 'acc_{dist1CrtDist}'),    ...
    ...
    getBetaMean(dynMdl{2}, 'acc_{dist1CrtDist}'),    ...
    ...
    getBetaMean(dynMdl{3}, 'acc_{dist1CrtDist}'),    ...
    ];

CI_col_dist_acc = [...
    getBetaCI(dynMdl{1}, 'acc_{dist1CrtDist}'),    ...
    ...
    getBetaCI(dynMdl{2}, 'acc_{dist1CrtDist}'),    ...
    ...
    getBetaCI(dynMdl{3}, 'acc_{dist1CrtDist}'),    ...
    ];


all_col_dist_acc{1} = getBetaAll(dynMdl{1}, 'acc_{dist1CrtDist}');
all_col_dist_acc{2} = getBetaAll(dynMdl{2}, 'acc_{dist1CrtDist}');
all_col_dist_acc{3} = getBetaAll(dynMdl{3}, 'acc_{dist1CrtDist}');


% n=nexttile([1,4]); hold on;

for ee = whichExps
    plot(repmat(x_col_dist_acc(ee), [npt(ee),1])+randn([npt(ee),1]).*.1.*ksdensity(all_col_dist_acc{ee},all_col_dist_acc{ee}), all_col_dist_acc{ee},   '.c');
end

errorbar(x_col_dist_acc(whichExps), M_col_dist_acc(whichExps), CI_col_dist_acc(whichExps),  'ok', 'MarkerFaceColor', 'c', 'LineWidth', 2, 'MarkerSize', 7);


% n.YLim = [0, Inf];
set(gca, 'TickDir', 'out', 'LineWidth', 1)
xlim([x_col_dist_acc(1)-.25, x_col_dist_acc(end) + .25])
yline(0, '-k', 'LineWidth', 1);

% yticks()
xticks([1, 1.75, 2.5])
xticklabels({'Exp1', 'Exp2', 'Exp3'})
title('dynamics betas')

ylim([-7.5, 6.5])
yticks([-7.5, 6.5])




% ============== saturate dynamics betas ==============
x_col_targ_acc = [1.5, 2.25];
x_col_dist_acc = [1, 1.75, 2.5]+.1;
x_col_sat_acc = [1, 1.5, 2.25]-.1;


M_col_dist_acc = [...
    getBetaMean(dynMdl{1}, 'acc_{dist1Crt}'),    ...
    ...
    getBetaMean(dynMdl{2}, 'acc_{dist1Crt}'),    ...
    ...
    getBetaMean(dynMdl{3}, 'acc_{dist1Crt}'),    ...
    ];

CI_col_dist_acc = [...
    getBetaCI(dynMdl{1}, 'acc_{dist1Crt}'),    ...
    ...
    getBetaCI(dynMdl{2}, 'acc_{dist1Crt}'),    ...
    ...
    getBetaCI(dynMdl{3}, 'acc_{dist1Crt}'),    ...
    ];


all_col_dist_acc{1} = getBetaAll(dynMdl{1}, 'acc_{dist1Crt}');
all_col_dist_acc{2} = getBetaAll(dynMdl{2}, 'acc_{dist1Crt}');
all_col_dist_acc{3} = getBetaAll(dynMdl{3}, 'acc_{dist1Crt}');


% n=nexttile([1,4]); hold on;

for ee = whichExps
    plot(repmat(x_col_sat_acc(ee), [npt(ee),1])+randn([npt(ee),1]).*.1.*ksdensity(all_col_dist_acc{ee},all_col_dist_acc{ee}), all_col_dist_acc{ee},   '.y');
end

errorbar(x_col_sat_acc(whichExps), M_col_dist_acc(whichExps), CI_col_dist_acc(whichExps),  'ok', 'MarkerFaceColor', 'y', 'LineWidth', 2, 'MarkerSize', 7);


% n.YLim = [0, Inf];
set(gca, 'TickDir', 'out', 'LineWidth', 1)
xlim([x_col_sat_acc(1)-.25, x_col_dist_acc(end) + .25])
yline(0, '-k', 'LineWidth', 1);

% yticks()
xticks([1, 1.75, 2.5])
xticklabels({'Exp1', 'Exp2', 'Exp3'})
title('dynamics betas')

ylim([-7.5, 6.5])
yticks([-7.5, 6.5])


%% compare color-motion adapt




disp_mdls = {'mot'}
for cc = 1:3
    
    for mm = 1:length(disp_mdls)
        
        % == color
        tbl = table;
        
        tbl.pname   = eval([disp_mdls{mm},'{cc}.results.paramName'])';
        tbl.df      = eval(['repmat(length(',disp_mdls{mm},'{cc}.data) - length(', disp_mdls{mm},'{cc}.results.stats.tval), [ length(', disp_mdls{mm},'{cc}.results.stats.tval), 1])']);
        tbl.t       = eval([disp_mdls{mm},'{cc}.results.stats.tval']);
        tbl.p       = eval([disp_mdls{mm},'{cc}.results.stats.p']);
        tbl.d       = eval([disp_mdls{mm},'{cc}.results.stats.alpha./sqrt(diag(', disp_mdls{mm}, '{cc}.results.stats.groupvar))']);
        
        fprintf('\n =============== exp %d %s =============== \n\n', cc, disp_mdls{mm});
        disp(tbl)
        
    end
    
    fprintf('\n\n\n\n\n');
    
end






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
    tcrt = ismember(col{cc}.results.paramName, 'acc_{distDist1}');
    mtcrt = ismember(mot{cc}.results.paramName, 'acc_{distDist1}');
    
    cv2 = col{cc}.results.stats.groupmeancovariance(tcrt,tcrt);
    cv3 = mot{cc}.results.stats.groupmeancovariance(mtcrt,mtcrt);
    
    accTarg_tval = (col{cc}.results.alpha(tcrt) - mot{cc}.results.alpha(mtcrt))./...
        sqrt(cv2 + cv3);
    
    accTarg_df = (col_idf*cv2 + mot_idf*cv3).^2 ./...
        (col_idf.*(col_idf*cv2).^2 + mot_idf.*(mot_idf*cv3).^2);
    
    accTarg_pval = 2*tcdf(abs(accTarg_tval), accTarg_df, 'upper');
    
    
    % acc targ
    tcrt = ismember(col{cc}.results.paramName, 'acc_{distDist1}');
    mtcrt = ismember(mot{cc}.results.paramName, 'acc_{targTarg1}');
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
    tcrt = ismember(col{cc}.results.paramName, 'rt_{distDist1}');
    mtcrt = ismember(mot{cc}.results.paramName, 'rt_{distDist1}');
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
    tcrt = ismember(col{cc}.results.paramName, 'rt_{distDist1}');
    mtcrt = ismember(mot{cc}.results.paramName, 'rt_{targTarg1}');
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
    fprintf('\n == mot-col distDist1 acc: t(%.5g)=%.5g, d=%.5g, p=%.5g',  accTarg_df, accTarg_tval, accTarg_tval./sqrt(accTarg_df), accTarg_pval)
    fprintf('\n == mot-col targTarg1 acc: t(%.5g)=%.5g, d=%.5g, p=%.5g',  accDist_df, accDist_tval, accDist_tval./sqrt(accDist_df), accDist_pval)
    fprintf('\n == mot-col distDist1 rt: t(%.5g)=%.5g, d=%.5g, p=%.5g',  rtTarg_df, rtTarg_tval, rtTarg_tval./sqrt(rtTarg_df), rtTarg_pval)
    fprintf('\n == mot-col targTarg1 rt: t(%.5g)=%.5g, d=%.5g, p=%.5g\n\n',  rtDist_df, rtDist_tval, rtDist_tval./sqrt(rtDist_df), rtDist_pval)
    
    if cc==1
        [accDist_tval, accDist_pval, accDist_df, rtDist_tval, rtDist_pval, rtDist_df] =deal(nan);
    end
    
    
    pvals(:,cc) = [accTarg_pval;  accDist_pval; rtTarg_pval; rtDist_pval];
    df(:,cc) = [accTarg_df; accDist_df; rtTarg_df; rtDist_df];
    
    sqrn(:,cc) = [cc_sqn; cc_sqn; cc_sqn; cc_sqn];
    sgnD(:,cc) = sign([accTarg_tval; accDist_tval; rtTarg_tval; rtDist_tval]);
    
    
end

sgnD = sgnD == sign(sgnD(:,end));


% combP = 2*normcdf(nansum(norminv(1-((realmin+pvals)./2)).*df, 2) ./ sqrt(nansum((df.^2),2)), 'upper')
combP = normcdf(nansum(sgnD.*-norminv(pvals).*sqn, 2) ./ sqrt(isfinite(pvals)*(sqn.^2)'), 'upper')



% df(cc)      = eval(['length(',disp_mdls{mm},'{cc}.data) - length(', disp_mdls{mm},'{cc}.results.stats.tval)']);
% sqn(cc) = sqrt(eval(['length(',disp_mdls{mm},'{cc}.data)']));
% 
% end
% 
% 
% dvals = squeeze(dstat(:,3,:));
% pvals = squeeze(dstat(:,2,:));
% sgnD = sign(dvals).*sign(dvals(:,end));
% 
% combP1 = 2*normcdf(nansum(norminv(1-((realmin+pvals)./2)).*df, 2) ./ sqrt(isfinite(pvals)*(df.^2)'), 'upper');
% combP = normcdf(nansum(sgnD.*-norminv(pvals).*sqn, 2) ./ sqrt(isfinite(pvals)*(sqn.^2)'), 'upper');






%% ====== REWARD ==========================================================
clear col polyMdl expMdl mot dynMdl


% tanh AR softplus
col{2} = load('2021-01-17_18-43_exp2_withinStatic.mat');
col{3} = load('2021-01-18_21-17_exp3_rewStatic.mat');
mot{3} = load('2021-03-06_05-07_exp3_rewMotionStatic.mat');

dynMdl{3} = load('2021-04-06_20-56_exp3_rewDualDyn.mat');


disp('');
disp('loaded');
disp('');



% plot stats
clc;
disp_mdls = {'col', 'mot', 'dynMdl'}
mdlNames = {'attend-color', 'attend-motion', 'dynamics'};

for cc = 3
    
    for mm = 1:length(disp_mdls)
        
        % transform parameters
        x = eval([disp_mdls{mm},'{cc}.results.stats.alpha']);
        for xx = 1:length(x)
           x(xx) =  eval([disp_mdls{mm},'{cc}.data(1).link{xx}(x(xx))']);
        end
        
        
        % correct parameter names
        n = eval([disp_mdls{mm},'{cc}.results.paramName'])';
        for nn = 1:length(n)
            n(nn) = strrep(n(nn), 'acc', 'choice');
            n(nn) = strrep(n(nn), 'rt_{choice}', 'rt_{acc}');
            n(nn) = strrep(n(nn), 'rt_{choiceTarg}', 'rt_{accTarg}');
            n(nn) = strrep(n(nn), 'rt_{choiceDist}', 'rt_{accDist}');
            n(nn) = strrep(n(nn), 'pwr', 'preGain');
            n(nn) = strrep(n(nn), 'dist1', 'prevdist');
            n(nn) = strrep(n(nn), 'Dist1', 'Prevdist');
            n(nn) = strrep(n(nn), 'targ1', 'prevtarg');
            n(nn) = strrep(n(nn), 'Targ1', 'Prevtarg');
            n(nn) = strrep(n(nn), 'crt', 'RT');
        end
        %         disp(n)
        
        % == make table
        tbl = table;
        
        tbl.paramName           = n;
        tbl.df                  = eval(['repmat(length(',disp_mdls{mm},'{cc}.data) - length(', disp_mdls{mm},'{cc}.results.stats.tval), [ length(', disp_mdls{mm},'{cc}.results.stats.tval), 1])']);
        tbl.paramVal            = x;
        tbl.tVal                = eval([disp_mdls{mm},'{cc}.results.stats.tval']);
        tbl.pVal                = eval([disp_mdls{mm},'{cc}.results.stats.p']);
        tbl.cohenD              = eval([disp_mdls{mm},'{cc}.results.stats.alpha./sqrt(diag(', disp_mdls{mm}, '{cc}.results.stats.groupvar))']);
        
        
        
        
        
        
        fprintf('\n =============== exp %d %s =============== \n\n', cc, disp_mdls{mm});
        %         disp(tbl)
        
        writetable(tbl, './tables/reward.xlsx', 'Sheet', sprintf('exp%d %s', cc, mdlNames{mm}))
        
    end
    
    %      fprintf('\n\n\n\n\n');
    
end





% ==== combined stats
aaallT = struct;
for mm = 1:length(disp_mdls)
    
    % combined table
    allpname = [];
    for cc = 3
        allpname =  [allpname; eval([disp_mdls{mm},'{cc}.results.paramName'])'];
    end
    unqpname = unique(upper(allpname));
    alluname{mm} = unqpname;
    
    [dstat] = nan(length(unqpname), 3, 3);
    
    clear df
    
    
    for cc = 3
        
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

    combP1 = 2*normcdf(nansum(norminv(1-((realmin+pvals)./2)).*df, 2) ./ sqrt(isfinite(pvals)*(df.^2)'), 'upper');
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





% ====== compare color & motion

for cc = 3
    
    npt = length(col{cc}.data);
    
    % acc targ
    tcrt = ismember(col{cc}.results.paramName, 'acc_{rewTarg}');
    cv2 = col{cc}.results.stats.groupmeancovariance(tcrt,tcrt);
    cv3 = mot{cc}.results.stats.groupmeancovariance(tcrt,tcrt);
    
    accTarg_tval = (col{cc}.results.alpha(tcrt) - mot{cc}.results.alpha(tcrt))./...
        sqrt(cv2 + cv3);
    
    accTarg_df = (cv2 + cv3).^2 ./...
        ((cv2.^2)./(npt-1) + (cv3.^2)./(npt-1));
    
    accTarg_pval = 2*tcdf(abs(accTarg_tval), accTarg_df, 'upper');
    
    
    % acc targ
    tcrt = ismember(col{cc}.results.paramName, 'acc_{rewDist}');
    cv2 = col{cc}.results.stats.groupmeancovariance(tcrt,tcrt);
    cv3 = mot{cc}.results.stats.groupmeancovariance(tcrt,tcrt);
    
    accDist_tval = (col{cc}.results.alpha(tcrt) - mot{cc}.results.alpha(tcrt))./...
        sqrt(cv2 + cv3);
    
    accDist_df = (cv2 + cv3).^2 ./...
        ((cv2.^2)./(npt-1) + (cv3.^2)./(npt-1));
    
    accDist_pval = 2*tcdf(abs(accDist_tval), accDist_df, 'upper');
    
    
    
    % rt targ
    tcrt = ismember(col{cc}.results.paramName, 'rt_{rewTarg}');
    cv2 = col{cc}.results.stats.groupmeancovariance(tcrt,tcrt);
    cv3 = mot{cc}.results.stats.groupmeancovariance(tcrt,tcrt);
    
    rtTarg_tval = (col{cc}.results.alpha(tcrt) - mot{cc}.results.alpha(tcrt))./...
        sqrt(cv2 + cv3);
    
    rtTarg_df = (cv2 + cv3).^2 ./...
        ((cv2.^2)./(npt-1) + (cv3.^2)./(npt-1));
    
    rtTarg_pval = 2*tcdf(abs(rtTarg_tval), rtTarg_df, 'upper');
    
    
    % rt targ
    tcrt = ismember(col{cc}.results.paramName, 'rt_{rewDist}');
    cv2 = col{cc}.results.stats.groupmeancovariance(tcrt,tcrt);
    cv3 = mot{cc}.results.stats.groupmeancovariance(tcrt,tcrt);
    
    rtDist_tval = (col{cc}.results.alpha(tcrt) - mot{cc}.results.alpha(tcrt))./...
        sqrt(cv2 + cv3);
    
    n = df(cc);
    df1 = ((cv2./n)+(cv3^2./n))^2;
    df2 = (cv2^2./(n^2*(n)))+(cv3^2./(n^2*(n)));
    rtDist_df = df1/df2;
    
%     rtDist_df = (cv2./df(cc) + cv3./df(cc)).^2 ./...
%         ((cv2./df(cc)).^2)./df(cc) + ((cv3./df(cc)).^2)./df(cc);
    
    rtDist_pval = 2*tcdf(abs(rtDist_tval), rtDist_df, 'upper');
    
    
    fprintf('\n exp %d', cc');
    fprintf('\n == col-mot rewTarg acc: t(%.5g)=%.5g, p=%.5g',  accTarg_df, accTarg_tval, accTarg_pval)
    fprintf('\n == col-mot rewDist acc: t(%.5g)=%.5g, p=%.5g',  accDist_df, accDist_tval, accDist_pval)
    fprintf('\n == col-mot rewTarg rt: t(%.5g)=%.5g, p=%.5g',  rtTarg_df, rtTarg_tval, rtTarg_pval)
    fprintf('\n == col-mot rewDist rt: t(%.5g)=%.5g, p=%.5g\n\n',  rtDist_df, rtDist_tval, rtDist_pval)
    
    
end












% ====== compare color & motion
clear pvals df rtDist_cvRatio sqrn sgnD
for cc = 3
    
    npt = length(col{cc}.data);
    cc_sqn = sqrt(npt);
    
    col_df  = npt - length(col{cc}.results.stats.tval);
    col_idf  = 1./col_df;
    mot_df  = npt - length(mot{cc}.results.stats.tval);
    mot_idf  = 1./mot_df;
    
    
    % acc targ
    tcrt = ismember(col{cc}.results.paramName, 'acc_{rewTarg}');
    mtcrt = ismember(mot{cc}.results.paramName, 'acc_{rewTarg}');
    
    cv2 = col{cc}.results.stats.groupmeancovariance(tcrt,tcrt);
    cv3 = mot{cc}.results.stats.groupmeancovariance(mtcrt,mtcrt);
    
    accTarg_tval = (col{cc}.results.alpha(tcrt) - mot{cc}.results.alpha(mtcrt))./...
        sqrt(cv2 + cv3);
    
    accTarg_df = (col_idf*cv2 + mot_idf*cv3).^2 ./...
        (col_idf.*(col_idf*cv2).^2 + mot_idf.*(mot_idf*cv3).^2);
    
    accTarg_pval = 2*tcdf(abs(accTarg_tval), accTarg_df, 'upper');
    
    
    % acc targ
    tcrt = ismember(col{cc}.results.paramName, 'acc_{rewDist}');
    mtcrt = ismember(mot{cc}.results.paramName, 'acc_{rewDist}');
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
    tcrt = ismember(col{cc}.results.paramName, 'rt_{rewTarg}');
    mtcrt = ismember(mot{cc}.results.paramName, 'rt_{rewTarg}');
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
    tcrt = ismember(col{cc}.results.paramName, 'rt_{rewDist}');
    mtcrt = ismember(mot{cc}.results.paramName, 'rt_{rewDist}');
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


%% plot reward






load('batlow.mat');
bat_cols = batlow;
bat_cIdx = round(linspace(1, length(bat_cols), 7));


mdl = dynMdl;


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

npt2 = [size(col{2}.data,2)];
npt = [size(col{3}.data,2)];

dyn_rt_tbl = array2table([rt, rt_acc, crt, rt_pred, rt_targ, rt_dist, rt_rew, rt_pt ], 'VariableNames', ...
    {'rt', 'acc', 'crt', 'rtY', 'targ', 'dist', 'rew', 'pt'});

dyn_acc_tbl = array2table([acc, cacc, acc_pred, acc_targ, acc_dist, acc_rt, acc_crt, d5rt, acc_rew, acc_pt], 'VariableNames', ...
    {'acc', 'cacc', 'accY', 'targ', 'dist', 'rt', 'crt', 'drt', 'rew', 'pt'});









mdl = col;


[rt, crt, rt_acc, rt_pred, rt_targ, rt_dist, rt_pt,...
    rt_rew, acc_rew,...
    acc, cacc, acc_pred, acc_targ, acc_dist, acc_pt,...
    all_ndt, exp2_ndt,...
    rt_blk, acc_blk, acc_sBlk, acc_rt, ...
    ] = deal([]);

pt_count = 1;

whichExps = 3;

nBlkDisc = 2;
nCohDisc = 11;


cols = [0,0,0; 255, 183, 3]/255;

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
        
        %                 acc_targ    = [acc_targ; discretize(tiedrank(out.data.targ0(out.acc.acc_sel)), nCohDisc)];
        %                 acc_dist    = [acc_dist; discretize(tiedrank((out.data.dist0(out.acc.acc_sel))), nCohDisc)];
        
        acc_targ    = [acc_targ; out.data.d11Targ0(out.acc.acc_sel)];
        acc_dist    = [acc_dist; -out.data.d11Dist0(out.acc.acc_sel)];
        
        if mm == 2
            acc_rew      = [acc_rew; out.data.d11Targ0(out.acc.acc_sel)*0];
        else
            acc_rew      = [acc_rew; grp2idx(out.data.rew(out.acc.acc_sel))];
        end
        
        acc_pt      = [acc_pt; out.acc.acc_X(:,1).*pt_count];
        
        
        pt_count = pt_count+1;
        
    end
    
    
end

npt2 = [size(col{2}.data,2)];
npt = [size(col{3}.data,2)];

rt_tbl = array2table([rt, rt_acc, crt, rt_pred, rt_targ, rt_dist, rt_rew, rt_pt ], 'VariableNames', ...
    {'rt', 'acc', 'crt', 'rtY', 'targ', 'dist', 'rew', 'pt'});

acc_tbl = array2table([acc, cacc, acc_pred, acc_targ, acc_dist, acc_rt, acc_rew, acc_pt], 'VariableNames', ...
    {'acc', 'cacc', 'accY', 'targ', 'dist', 'rt', 'rew', 'pt'});















% plot
figure('units','inch','position',[0,0,10,8]); hold on;
tiledlayout(2,20, 'TileSpacing', 'compact', 'Padding', 'compact');





% =================== TARGET

dotOff = linspace(-.15,.15,nBlkDisc);


% ============== targ acc ==============
grpTbl = grpstats(acc_tbl, {'pt','targ', 'rew'});
grpTbl = grpstats(grpTbl, {'targ', 'rew'}, {'mean','sem'});

nexttile([1,8]); hold on;



blkSel = grpTbl.rew == 0;
% plot([1:nCohDisc], grpTbl.mean_mean_acc(blkSel),...
%     'ok','MarkerFaceColor', 'w','LineWidth', 1, 'MarkerSize', 5);

allY = [];
for ii = 1:nBlkDisc
    
    blkSel = grpTbl.rew == ii;
    
    e=errorbar([1:nCohDisc] + dotOff(ii), grpTbl.mean_mean_acc(blkSel), grpTbl.sem_mean_cacc(blkSel),...
        'o', 'Color', cols(ii,:), 'MarkerFaceColor', cols(ii,:), 'LineWidth', 1);
    
    plot([1:nCohDisc] + dotOff(ii), grpTbl.mean_mean_accY(blkSel), '-', 'LineWidth', 2, 'Color', cols(ii,:));
    
    allY = [allY; e.YData(:)];
    
end

yticks(ylim)

title('Target \times Reward (Accuracy)')
set(gca, 'TickDir', 'out', 'LineWidth', 1)
xticks([])
xlabel('Target Coherence')



% ============== targ rt ==============
grpTbl = grpstats(rt_tbl, {'pt','targ', 'rew'}, {'mean','median'});
grpTbl = grpstats(grpTbl, {'targ', 'rew'}, {'mean','sem'});

nexttile([1,8]); hold on;


blkSel = grpTbl.rew == 0;
%  plot([1:nCohDisc], exp(grpTbl.mean_mean_rt(blkSel)) + mean(exp2_ndt),...
%         'ok','MarkerFaceColor', 'w','LineWidth', 1, 'MarkerSize', 5);
allY =[];
for ii = 1:nBlkDisc
    
    blkSel = grpTbl.rew == ii;
    
    e=errorbar([1:nCohDisc] + dotOff(ii), exp(grpTbl.mean_mean_rt(blkSel)) + mean(all_ndt), exp(grpTbl.mean_mean_rt(blkSel)).*(grpTbl.sem_mean_crt(blkSel)),...
        'o', 'Color', cols(ii,:), 'MarkerFaceColor', cols(ii,:), 'LineWidth', 1);
    
    plot([1:nCohDisc] + dotOff(ii), exp(grpTbl.mean_mean_rtY(blkSel)) + mean(all_ndt), '-', 'LineWidth', 2, 'Color', cols(ii,:));
    allY = [allY; e.YData(:)];
    
end

yticks(ylim)

title('Target \times Reward (RT)')
set(gca, 'TickDir', 'out', 'LineWidth', 1)
xticks([])
xlabel('Target Coherence')









% ============== targ betas ==============
x_col_targ = [1, 1.5];


M_col_dist_acc = [...
    getBetaMean(col{3}, 'acc_{rewTarg}'),    ...
    ];

CI_col_dist_acc = [...
    getBetaCI(col{3}, 'acc_{rewTarg}'),    ...
    ];


all_col_dist_acc = getBetaAll(col{3}, 'acc_{rewTarg}');



M_col_dist_rt = [...
    getBetaMean(col{3}, 'rt_{rewTarg}'),    ...
    ];

CI_col_dist_rt = [...
    getBetaCI(col{3}, 'rt_{rewTarg}'),    ...
    ];


all_col_dist_rt = getBetaAll(col{3}, 'rt_{rewTarg}');





n=nexttile([1,4]); hold on;

yyaxis('left');
plot(repmat(x_col_targ(1), [npt,1])+randn([npt,1]).*.0075.*ksdensity(all_col_dist_acc,all_col_dist_acc), all_col_dist_acc,   '.', 'color', rtacc_cols{1});
errorbar(x_col_targ(1), M_col_dist_acc, CI_col_dist_acc,  'ok', 'MarkerFaceColor', 'k', 'LineWidth', 2, 'MarkerSize', 7);

ylim([-.4, 1])
yticks(ylim)


yyaxis('right');
plot(repmat(x_col_targ(2), [npt,1])+randn([npt,1]).*.005.*ksdensity(all_col_dist_rt,all_col_dist_rt), all_col_dist_rt,   '.', 'color', rtacc_cols{2});
errorbar(x_col_targ(2), M_col_dist_rt, CI_col_dist_rt,  'ok', 'MarkerFaceColor', 'k', 'LineWidth', 2, 'MarkerSize', 7);

ylim([-.2, 0.5])
yticks(ylim)


% n.YLim = [0, Inf];
set(gca, 'TickDir', 'out', 'LineWidth', 1)
xlim([x_col_targ(1)-.25, x_col_targ(end) + .25])
yline(0, '-k', 'LineWidth', 1);
xticks([])

title('target betas')












% =================== DISTRACTOR


% ============== distarctor acc ==============
grpTbl = grpstats(acc_tbl, {'pt','dist', 'rew'});
grpTbl = grpstats(grpTbl, {'dist', 'rew'}, {'mean','sem'});

nexttile([1,8]); hold on;

blkSel = grpTbl.rew == 0;
%  plot([1:nCohDisc], grpTbl.mean_mean_acc(blkSel),...
%         'ok','MarkerFaceColor', 'w','LineWidth', 1, 'MarkerSize', 5);
allY = [];
for ii = 1:nBlkDisc
    
    blkSel = grpTbl.rew == ii;
    
    e=errorbar([1:nCohDisc] + dotOff(ii), grpTbl.mean_mean_acc(blkSel), grpTbl.sem_mean_cacc(blkSel),...
        'o', 'Color', cols(ii,:), 'MarkerFaceColor', cols(ii,:), 'LineWidth', 1);
    
    plot([1:nCohDisc] + dotOff(ii), grpTbl.mean_mean_accY(blkSel), '-', 'LineWidth', 2, 'Color', cols(ii,:));
    allY = [allY; e.YData(:)];
    
end

yticks(ylim)

title('Distractor \times Reward (Accuracy)')
set(gca, 'TickDir', 'out', 'LineWidth', 1)
xticks([])
xlabel('Distractor Congruence')



% ============== distarctor rt ==============
grpTbl = grpstats(rt_tbl, {'pt','dist', 'rew'}, {'mean','median'});
grpTbl = grpstats(grpTbl, {'dist', 'rew'}, {'mean','sem'});

nexttile([1,8]); hold on;

blkSel = grpTbl.rew == 0;
%  plot([1:nCohDisc], exp(grpTbl.mean_mean_rt(blkSel)) + mean(exp2_ndt),...
%         'ok','MarkerFaceColor', 'w','LineWidth', 1, 'MarkerSize', 5);
allY=[];
for ii = 1:nBlkDisc
    
    blkSel = grpTbl.rew == ii;
    
    e=errorbar([1:nCohDisc] + dotOff(ii), exp(grpTbl.mean_mean_rt(blkSel)) + mean(all_ndt), exp(grpTbl.mean_mean_rt(blkSel)).*(grpTbl.sem_mean_crt(blkSel)),...
        'o', 'Color', cols(ii,:), 'MarkerFaceColor', cols(ii,:), 'LineWidth', 1);
    
    plot([1:nCohDisc] + dotOff(ii), exp(grpTbl.mean_mean_rtY(blkSel)) + mean(all_ndt), '-', 'LineWidth', 2, 'Color', cols(ii,:));
    
    allY = [allY; e.YData(:)];
    
end

yticks(ylim)

title('Distractor \times Reward (RT)')
set(gca, 'TickDir', 'out', 'LineWidth', 1)
xticks([])
xlabel('Distractor Congruence')




% ============== targ dynamics betas ==============
% ============== targ betas ==============
x_col_targ = [1, 1.5];


M_col_dist_acc = [...
    getBetaMean(col{3}, 'acc_{rewDist}'),    ...
    ];

CI_col_dist_acc = [...
    getBetaCI(col{3}, 'acc_{rewDist}'),    ...
    ];


all_col_dist_acc = getBetaAll(col{3}, 'acc_{rewDist}');



M_col_dist_rt = [...
    getBetaMean(col{3}, 'rt_{rewDist}'),    ...
    ];

CI_col_dist_rt = [...
    getBetaCI(col{3}, 'rt_{rewDist}'),    ...
    ];


all_col_dist_rt = getBetaAll(col{3}, 'rt_{rewDist}');





n=nexttile([1,4]); hold on;

yyaxis('left');
plot(repmat(x_col_targ(1), [npt,1])+randn([npt,1]).*.005.*ksdensity(all_col_dist_acc,all_col_dist_acc), all_col_dist_acc,   '.', 'color', rtacc_cols{1});
errorbar(x_col_targ(1), M_col_dist_acc, CI_col_dist_acc,  'ok', 'MarkerFaceColor', 'k', 'LineWidth', 2, 'MarkerSize', 7);

ylim([-.2, .35])
yticks(ylim)


yyaxis('right');
plot(repmat(x_col_targ(2), [npt,1])+randn([npt,1]).*.0005.*ksdensity(all_col_dist_rt,all_col_dist_rt), all_col_dist_rt,   '.', 'color', rtacc_cols{2});
errorbar(x_col_targ(2), M_col_dist_rt, CI_col_dist_rt,  'ok', 'MarkerFaceColor', 'k', 'LineWidth', 2, 'MarkerSize', 7);

ylim([-.02, 0.035])
yticks(ylim)


% n.YLim = [0, Inf];
set(gca, 'TickDir', 'out', 'LineWidth', 1)
xlim([x_col_targ(1)-.25, x_col_targ(end) + .25])
yline(0, '-k', 'LineWidth', 1);
xticks([])









figure('units','inch','position',[0,0,1.5*7.5,1.5*2]); hold on;
tiledlayout(1,5, 'TileSpacing', 'none', 'Padding', 'none');


% ============== plot dynamics effects ==============



grpTbl = grpstats(dyn_acc_tbl, {'pt', 'targ', 'drt' 'rew'});
grpTbl = grpstats(grpTbl, {'targ', 'drt', 'rew'}, {'mean','sem'});


nexttile([1,2]); hold on;
yline(0, '--k', 'LineWidth', 1)

xxoff = linspace(-.1,.1,5);
for ii = [1:5]
    
    plot([1:5]+xxoff(ii), grpTbl.mean_mean_acc(grpTbl.drt == ii & grpTbl.rew == 2)-grpTbl.mean_mean_acc(grpTbl.drt == ii & grpTbl.rew == 1),...
        'o', 'LineWidth', 1, 'color', bat_cols(bat_cIdx(ii+1),:),  'MarkerFaceColor', bat_cols(bat_cIdx(ii+1),:));
    
       pltp(ii) = plot([1:5]+xxoff(ii), grpTbl.mean_mean_accY(grpTbl.drt == ii & grpTbl.rew == 2)-grpTbl.mean_mean_accY(grpTbl.drt == ii & grpTbl.rew == 1),...
        '-', 'LineWidth', 2, 'color', bat_cols(bat_cIdx(ii+1),:));
    
%         accii = grpTbl.mean_mean_acc(grpTbl.drt == ii & grpTbl.rew == 2)-grpTbl.mean_mean_acc(grpTbl.drt == ii & grpTbl.rew == 1);
%         accyii = grpTbl.mean_mean_accY(grpTbl.drt == ii & grpTbl.rew == 2)-grpTbl.mean_mean_accY(grpTbl.drt == ii & grpTbl.rew == 1);
%         for jj = 1:length(accii)
%              plot([jj,jj]+xxoff(ii), [accii(jj), accyii(jj)],...
%             '-', 'LineWidth', .25, 'color', bat_cols(bat_cIdx(ii+1),:));
%         end
          
end

ylim([-.02, .09])
yticks(ylim)
xlim([.5,5.5])
title('reward interaction')
ylabel('HI - NO rew (accuracy)')
xlabel('target coherence')
set(gca, 'TickDir', 'out', 'LineWidth', 1);
% legend(pltp, {'fast','','med','','slow'}, 'Location', 'northeast')

title('dist betas')














% ======  plot dynamics ==================================================
nptSel = [ones(npt2,1)*2; ones(npt,1)*3];
grpTbl = grpstats(acc_tbl, {'pt'});
rtMx = repmat([.5:.001:1], [size(grpTbl,1),1]) - grpTbl.mean_rt;

mRT   = mean(acc_tbl.rt);
rtRange = [.5:.001:1] - mRT;




% plot
n=nexttile([1,2]); hold on;

xline(prctile(acc_rt(acc_rew~=0), 50), '-k');
xline(prctile(acc_rt(acc_rew~=0), 25), '--k');
xline(prctile(acc_rt(acc_rew~=0), 75), '--k');

rtRange = [.5:.001:1];


darkCol = linspace(1,.2,nBlkDisc);

rews = [-1 1];
conRew = acc_rew;
conRew(conRew==1) = -1;
conRew(conRew==2) = 1;


for bb = [1:nBlkDisc]
    
    mBlk = rews(bb);
    
    clear mRT
    for pp = 1:size(grpTbl,1)
        
        if std(conRew(acc_pt==pp)) ~= 0
%             mRT(pp,1) = (conRew(acc_pt==pp)\center(acc_tbl.rt(acc_pt==pp)))*mBlk + mean(acc_tbl.rt(acc_pt==pp));
            mRT(pp,1) = nanmean(dyn_acc_tbl.crt(dyn_acc_tbl.pt== (pp) & dyn_acc_tbl.rew == bb)) + mean(acc_tbl.rt(acc_pt==pp));
        else
        end
        
    end
    rtMx = repmat([.5:.001:1], [size(grpTbl,1),1]) - mRT;
    
    
    
    
    
    dyn_lapse  = [...
        1./(1+exp(-(getBetaAll(dynMdl{3}, 'acc_{lapse}') + rtMx.*getBetaAll(dynMdl{3}, 'acc_{crt}') + mBlk.*getBetaAll(dynMdl{3}, 'acc_{rew}') + mBlk.*rtMx.*getBetaAll(dynMdl{3}, 'acc_{rewCrt}'))));
        ];
    
    
    
    dyn_targ  = (1-dyn_lapse).*[...
        getBetaAll(dynMdl{3}, 'acc_{targ}') + rtMx.*getBetaAll(dynMdl{3}, 'acc_{crtTarg}') + mBlk.*getBetaAll(dynMdl{3}, 'acc_{rewTarg}') + mBlk.*rtMx.*getBetaAll(dynMdl{3}, 'acc_{rewCrtTarg}');...
        ];
    
    
    dyn_dist = (1-dyn_lapse).*[...
        getBetaAll(dynMdl{3}, 'acc_{dist}') + rtMx.*getBetaAll(dynMdl{3}, 'acc_{crtDist}') + mBlk.*getBetaAll(dynMdl{3}, 'acc_{rewDist}') + mBlk.*rtMx.*getBetaAll(dynMdl{3}, 'acc_{rewCrtDist}');...
        ];
    
    
    
    
    targN = size(dyn_targ,1);
    dynN = size(dyn_dist,1);
    
    
    
    
    
    yline(0, '--k', 'lineWidth', 3);
    
    
    %     pltT=plot(rtRange, mean((dyn_targ)), '-', 'LineWidth', 3, 'Color', cols(cIdx(bb+1),:));
    %     plot(rtRange, mean((dyn_targ)) + (std((dyn_targ))./sqrt(targN)), '-', 'Color', cols(cIdx(bb+1),:), 'LineWidth', .5)
    %     plot(rtRange, mean((dyn_targ)) - (std((dyn_targ))./sqrt(targN)), '-', 'Color', cols(cIdx(bb+1),:), 'LineWidth', .5)
    
    pltT=plot(rtRange, mean((dyn_targ)), '-', 'LineWidth', 3, 'Color', [0, darkCol(bb), 0]);
    plot(rtRange, mean((dyn_targ)) + (std((dyn_targ))./sqrt(targN)), '-', 'Color', [0, darkCol(bb), 0], 'LineWidth', .5)
    plot(rtRange, mean((dyn_targ)) - (std((dyn_targ))./sqrt(targN)), '-', 'Color', [0, darkCol(bb), 0], 'LineWidth', .5)
    
    
    %     pltD=plot(rtRange, mean((dyn_dist)), '-', 'LineWidth', 3, 'Color', cols(cIdx(bb+1),:));
    %     plot(rtRange, mean((dyn_dist)) + std((dyn_dist))./sqrt(dynN), '-', 'Color', cols(cIdx(bb+1),:),  'LineWidth', .5)
    %     plot(rtRange, mean((dyn_dist)) - std((dyn_dist))./sqrt(dynN), '-', 'Color', cols(cIdx(bb+1),:), 'LineWidth', .5)
    
    
    pltD=plot(rtRange, mean((dyn_dist)), '-', 'LineWidth', 3, 'Color', [0, darkCol(bb), darkCol(bb)]);
    plot(rtRange, mean((dyn_dist)) + std((dyn_dist))./sqrt(dynN), '-', 'Color', [0, darkCol(bb), darkCol(bb)],  'LineWidth', .5)
    plot(rtRange, mean((dyn_dist)) - std((dyn_dist))./sqrt(dynN), '-', 'Color', [0, darkCol(bb), darkCol(bb)], 'LineWidth', .5)
    
    
    
end


% xlim([min(rtRange+mRT), max(rtRange+mRT)])
% ylim([-8, 6])
% yticks([-8, 6])

xlabel('time (s)')
ylabel('coherence sensitivity')
title('sensitivity dynamics')
set(gca, 'TickDir', 'out', 'LineWidth', 1)

% legend([pltT, pltD], {'target', 'distractor'}, 'Location', 'best')













% ============== targ dynamics betas ==============
x_col_targ_acc = [1.5, 2.25];
x_col_dist_acc = [1, 1.75, 2.5]-.1;


M_col_dist_acc = [...
    %     getBetaMean(polyMdl{1}, 'acc_{blkCrtTarg}'),    ...
    
    ...
    getBetaMean(dynMdl{3}, 'acc_{rewCrtTarg}'),    ...
    ];

CI_col_dist_acc = [...
    ...
    getBetaCI(dynMdl{3}, 'acc_{rewCrtTarg}'),    ...
    ];


% all_col_dist_acc{1} = getBetaAll(polyMdl{1}, 'acc_{blkCrtTarg}');
all_col_dist_acc = getBetaAll(dynMdl{3}, 'acc_{rewCrtTarg}');


n=nexttile([1,1]); hold on;

plot(repmat(x_col_dist_acc(2), [npt(1),1])+randn([npt(1),1]).*.1.*ksdensity(all_col_dist_acc,all_col_dist_acc), all_col_dist_acc,   '.g');

errorbar(x_col_dist_acc(2), M_col_dist_acc(1), CI_col_dist_acc(1),  'ok', 'MarkerFaceColor', 'k', 'LineWidth', 2, 'MarkerSize', 7);


% n.YLim = [0, Inf];
set(gca, 'TickDir', 'out', 'LineWidth', 1)
xlim([x_col_dist_acc(2)-.25, x_col_dist_acc(2) + .25])
yline(0, '-k', 'LineWidth', 1);
% ylim([-8, 6])
% yticks([-8, 6])



% ============== dist dynamics betas ==============
x_col_targ_acc = [1.5, 2.25];
x_col_dist_acc = [1, 1.75, 2.5]+.1;


M_col_dist_acc = [...
    ...
    getBetaMean(dynMdl{3}, 'acc_{rewCrtDist}'),    ...
    ];

CI_col_dist_acc = [...
    ...
    getBetaCI(dynMdl{3}, 'acc_{rewCrtDist}'),    ...
    ];


all_col_dist_acc = getBetaAll(dynMdl{3}, 'acc_{rewCrtDist}');


% n=nexttile([1,4]); hold on;

plot(repmat(x_col_dist_acc(2), [npt(1),1])+randn([npt(1),1]).*.1.*ksdensity(all_col_dist_acc,all_col_dist_acc), all_col_dist_acc,   '.c');
errorbar(x_col_dist_acc(2), M_col_dist_acc, CI_col_dist_acc,  'ok', 'MarkerFaceColor', 'k', 'LineWidth', 2, 'MarkerSize', 7);


% n.YLim = [0, Inf];
set(gca, 'TickDir', 'out', 'LineWidth', 1)
xlim([x_col_dist_acc(2)-.1-.25, x_col_dist_acc(2)-.1 + .25])
yline(0, '-k', 'LineWidth', 1);


xticks([1.75, 2.5])
xticklabels({'Exp3'})

% ylim([-8, 6])
% yticks([-8, 6])








% ============== saturate dynamics betas ==============
x_col_targ_acc = [1.5, 2.25];
x_col_dist_acc = [1, 1.75, 2.5]+.1;
x_col_sat_acc = [1, 1.5, 2.25]-.1;


M_col_sat_acc = [...
    getBetaMean(dynMdl{3}, 'acc_{rewCrt}'),    ...
    ];

CI_col_sat_acc = [...
    getBetaCI(dynMdl{3}, 'acc_{rewCrt}'),    ...
    ];


all_col_sat_acc = getBetaAll(dynMdl{3}, 'acc_{rewCrt}');




plot(repmat(x_col_sat_acc(2), [npt(1),1])+randn([npt(1),1]).*.1.*ksdensity(all_col_sat_acc,all_col_sat_acc), all_col_sat_acc,   '.y');
errorbar(x_col_sat_acc(2), M_col_sat_acc, CI_col_sat_acc,  'ok', 'MarkerFaceColor', 'k', 'LineWidth', 2, 'MarkerSize', 7);


% n.YLim = [0, Inf];
set(gca, 'TickDir', 'out', 'LineWidth', 1)
xlim([x_col_sat_acc(2)-.25, x_col_dist_acc(2) + .25])
yline(0, '-k', 'LineWidth', 1);

% yticks()
xticks([1, 1.75, 2.5])
xticklabels({'Exp1', 'Exp2', 'Exp3'})
title('dynamics betas')

ylim([-6.1, 5])
yticks([-6, 5])









%%  %% ====== LOAD DOUBLE DISSOCIATION BT REWARD & ADAPTATION ===================
clear rew adapt



% ==== LOAD

% tanh AR softplus
rew = load('2021-01-18_21-17_exp3_rewStatic.mat');


% tanh AR softplus
adapt{1} = load('2021-01-18_02-26_exp1_adaptStatic.mat');
adapt{2} = load('2021-01-20_05-44_exp2_adaptStatic.mat');
adapt{3} = load('2021-01-18_22-53_exp3_adaptStatic.mat');


%% ==== RUN DOUBLE DISSOCIATION BT REWARD & ADAPTATION ===================

npt2 = length(adapt{2}.data);
snpt2 = sqrt(npt2);
a2df  = npt2 - length(adapt{2}.results.stats.tval);


npt3 = length(adapt{3}.data);
snpt3 = sqrt(npt3);
a3df  = npt3 - length(adapt{3}.results.stats.tval);




% idx
rewTargIdx  = ismember(rew.results.paramName, 'acc_{rewTarg}');
rewDistIdx  = ismember(rew.results.paramName, 'acc_{rewDist}');
adaptTargIdx  = ismember(adapt{2}.results.paramName, 'acc_{targDist1}');
adaptDistIdx  = ismember(adapt{2}.results.paramName, 'acc_{distDist1}');

% mean
rewTargMu = rew.results.alpha(rewTargIdx);
rewDistMu = rew.results.alpha(rewDistIdx);
a2TargMu = adapt{2}.results.alpha(adaptTargIdx);
a2DistMu = adapt{2}.results.alpha(adaptDistIdx);
a3TargMu = adapt{3}.results.alpha(adaptTargIdx);
a3DistMu = adapt{3}.results.alpha(adaptDistIdx);

mu = abs([rewTargMu, rewDistMu, a2TargMu, a2DistMu, a3TargMu, a3DistMu]);

% cov
cv = blkdiag(...
    rew.results.stats.groupmeancovariance(rewTargIdx|rewDistIdx, rewTargIdx|rewDistIdx),...
    adapt{2}.results.stats.groupmeancovariance(adaptTargIdx|adaptDistIdx, adaptTargIdx|adaptDistIdx),...
    adapt{3}.results.stats.groupmeancovariance(adaptTargIdx|adaptDistIdx, adaptTargIdx|adaptDistIdx));

gv = blkdiag(...
    rew.results.stats.groupvar(rewTargIdx|rewDistIdx, rewTargIdx|rewDistIdx),...
    adapt{2}.results.stats.groupvar(adaptTargIdx|adaptDistIdx, adaptTargIdx|adaptDistIdx),...
    adapt{3}.results.stats.groupvar(adaptTargIdx|adaptDistIdx, adaptTargIdx|adaptDistIdx));

scv = (diag(cv))';
sgv = sqrt(diag(gv))';

% contrasts

% === within-reward
desc = {'within-rew', 'within-adapt', 'within-target','within-distractor', 'interaction'}
cons = [...
    [snpt2+snpt3, -(snpt2+snpt3), 0, 0, 0, 0]./(snpt2+snpt3);
    [0, 0, -snpt2, snpt2, -snpt3, snpt3]./(snpt2+snpt3);
    [snpt2+snpt3, 0, -snpt2, 0, -snpt3, 0]./(snpt2+snpt3);
    [0, -(snpt2+snpt3), 0, snpt2, 0, snpt3]./(snpt2+snpt3);
    [snpt2+snpt3, -(snpt2+snpt3), -snpt2, snpt2, -snpt3, snpt3]./(snpt2+snpt3)...
    ]


for ii = 1:length(desc)
    
    fprintf('\n\n%s\n', desc{ii})
    
    tval = (cons(ii,:)*mu') ./ sqrt(cons(ii,:)*cv*cons(ii,:)')
    
    vsel = cons(ii,:)~=0;
    v   = [a3df, a3df, a2df, a2df, a3df, a3df];
    iv  = 1./v;
    
    df = sum(iv(vsel).*scv(vsel)).^2 ./...
        sum(((iv(vsel).*scv(vsel)).^2)./v(vsel))
    
    pval = 2*tcdf(abs(tval), df, 'upper')
    
    
end







% plot
figure('units','inch','position',[0,0,16,8]); hold on;
tiledlayout(1,2, 'TileSpacing', 'compact', 'Padding', 'compact');



% errorbar([.9 1.4], mu(3:4), sqrt(scv(3:4)), '-or', 'MarkerFaceColor', 'r', 'LineWidth', 1, 'MarkerSize', 9);
% papt=errorbar([1,1.5], mu(5:6), sqrt(scv(5:6)), '-or', 'MarkerFaceColor', 'r', 'LineWidth', 1, 'MarkerSize', 9);
% prew=errorbar([1.1,1.6], mu(1:2), sqrt(scv(1:2)), '-ob', 'MarkerFaceColor', 'b', 'LineWidth', 1, 'MarkerSize', 9);

nexttile; hold on;
errorbar([1,1.5]-.01, mu(3:4), sqrt(scv(3:4)), '-or', 'MarkerFaceColor', 'r', 'LineWidth', 1.25, 'MarkerSize', 9);
papt=errorbar([1,1.5], mu(5:6), sqrt(scv(5:6)), '-or', 'MarkerFaceColor', 'r', 'LineWidth', 1.25, 'MarkerSize', 9);
prew=errorbar([1,1.5]+.01, mu(1:2), sqrt(scv(1:2)), '-ob', 'MarkerFaceColor', 'b', 'LineWidth', 1.25, 'MarkerSize', 9);


yline(0, '--k');

legend([prew,papt], {'reward', 'adaptation'}, 'Location', 'NorthWest')
xlim([.8, 1.7])
xticks([1,1.5])
xticklabels({'target', 'distractor'})
ylabel('absolute regression weight')
set(gca, 'TickDir', 'out', 'LineWidth', 1)
title('adaptation-reward dissociation  (Accuracy)')
ylim([-.1, .6])







% idx
rewTargIdx  = ismember(rew.results.paramName, 'rt_{rewTarg}');
rewDistIdx  = ismember(rew.results.paramName, 'rt_{rewDist}');
adaptTargIdx  = ismember(adapt{2}.results.paramName, 'rt_{targDist1}');
adaptDistIdx  = ismember(adapt{2}.results.paramName, 'rt_{distDist1}');

% mean
rewTargMu = rew.results.alpha(rewTargIdx);
rewDistMu = rew.results.alpha(rewDistIdx);
a2TargMu = adapt{2}.results.alpha(adaptTargIdx);
a2DistMu = adapt{2}.results.alpha(adaptDistIdx);
a3TargMu = adapt{3}.results.alpha(adaptTargIdx);
a3DistMu = adapt{3}.results.alpha(adaptDistIdx);

mu = abs([rewTargMu, rewDistMu, a2TargMu, a2DistMu, a3TargMu, a3DistMu]);

% cov
cv = blkdiag(...
    rew.results.stats.groupmeancovariance(rewTargIdx|rewDistIdx, rewTargIdx|rewDistIdx),...
    adapt{2}.results.stats.groupmeancovariance(adaptTargIdx|adaptDistIdx, adaptTargIdx|adaptDistIdx),...
    adapt{3}.results.stats.groupmeancovariance(adaptTargIdx|adaptDistIdx, adaptTargIdx|adaptDistIdx));

gv = blkdiag(...
    rew.results.stats.groupvar(rewTargIdx|rewDistIdx, rewTargIdx|rewDistIdx),...
    adapt{2}.results.stats.groupvar(adaptTargIdx|adaptDistIdx, adaptTargIdx|adaptDistIdx),...
    adapt{3}.results.stats.groupvar(adaptTargIdx|adaptDistIdx, adaptTargIdx|adaptDistIdx));

scv = (diag(cv))';
sgv = sqrt(diag(gv))';




% === within-reward
desc = {'within-rew', 'within-adapt', 'within-target','within-distractor', 'interaction'}
cons = [...
    [snpt2+snpt3, -(snpt2+snpt3), 0, 0, 0, 0]./(snpt2+snpt3);
    [0, 0, -snpt2, snpt2, -snpt3, snpt3]./(snpt2+snpt3);
    [snpt2+snpt3, 0, -snpt2, 0, -snpt3, 0]./(snpt2+snpt3);
    [0, -(snpt2+snpt3), 0, snpt2, 0, snpt3]./(snpt2+snpt3);
    [snpt2+snpt3, -(snpt2+snpt3), -snpt2, snpt2, -snpt3, snpt3]./(snpt2+snpt3)...
    ]


for ii = 1:length(desc)
    
    fprintf('\n\n%s\n', desc{ii})
    
    tval = (cons(ii,:)*mu') ./ sqrt(cons(ii,:)*cv*cons(ii,:)')
    
    vsel = cons(ii,:)~=0;
    v   = [a3df, a3df, a2df, a2df, a3df, a3df];
    iv  = 1./v;
    
    df = sum(iv(vsel).*scv(vsel)).^2 ./...
        sum(((iv(vsel).*scv(vsel)).^2)./v(vsel))
    
    pval = 2*tcdf(abs(tval), df, 'upper')
    
    
end










nexttile; hold on;
% errorbar([.9 1.4], mu(3:4), sqrt(scv(3:4)), '-or', 'MarkerFaceColor', 'r', 'LineWidth', 1, 'MarkerSize', 9);
% papt=errorbar([1,1.5], mu(5:6), sqrt(scv(5:6)), '-or', 'MarkerFaceColor', 'r', 'LineWidth', 1, 'MarkerSize', 9);
% prew=errorbar([1.1,1.6], mu(1:2), sqrt(scv(1:2)), '-ob', 'MarkerFaceColor', 'b', 'LineWidth', 1, 'MarkerSize', 9);
errorbar([1 1.5]-.01, mu(3:4), sqrt(scv(3:4)), '-or', 'MarkerFaceColor', 'r', 'LineWidth', 1.25, 'MarkerSize', 9);
papt=errorbar([1,1.5], mu(5:6), sqrt(scv(5:6)), '-or', 'MarkerFaceColor', 'r', 'LineWidth', 1.25, 'MarkerSize', 9);
prew=errorbar([1,1.5]+.01, mu(1:2), sqrt(scv(1:2)), '-ob', 'MarkerFaceColor', 'b', 'LineWidth', 1.25, 'MarkerSize', 9);
yline(0, '--k');

legend([prew,papt], {'reward', 'adaptation'}, 'Location', 'NorthWest')
xlim([.8, 1.7])
xticks([1,1.5])
xticklabels({'target', 'distractor'})
ylabel('absolute regression weight')
set(gca, 'TickDir', 'out', 'LineWidth', 1)
title('adaptation-reward dissociation (RT)')
ylim([-.01, .06])










