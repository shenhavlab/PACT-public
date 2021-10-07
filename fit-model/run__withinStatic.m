%% run__est_withinStatic =============================================
%
% function run__est_withinStatic(exp)
%
%
%
%
% Fit RDM non-linear regression using EM
%
%
%  harrison ritz (harrison.ritz@gmail.com) 2020
%
% =========================================================================


function run__withinStatic(name)

serial = 0;
if nargin < 1
    
    serial = 1;
    
    name = 'exp2';
    
end


saveAs = sprintf('%s_withinStatic', name)
likfun = 'lik_withinStatic'


%% add paths

addpath(genpath('../util'));


%% load data

data = load_RDM(name);



% select trials

minRT       = 0.20
maxRT       = 2.00

for pp = 1:length(data)
    
    
    data(pp).rt_sel = ...
        (data(pp).acc0 == 1) & ...
        (data(pp).acc1 == 1) & ...
        (data(pp).rt0 > minRT) & ...
        (data(pp).rt0 < maxRT) & ...
        (data(pp).rt1 > minRT) & ...
        (data(pp).rt1 < maxRT) & ...
        ismember(data(pp).ordType, 'rand') & ...
        ismember(data(pp).colBlk, 'color');
    
    
    data(pp).acc_sel = ...
        (data(pp).acc1 == 1) & ...
        (data(pp).rt0 > minRT) & ...
        (data(pp).rt0 < maxRT) & ...
        (data(pp).rt1 > minRT) & ...
        (data(pp).rt1 < maxRT) & ...
        ismember(data(pp).ordType, 'rand') & ...
        ismember(data(pp).colBlk, 'color');
    
    
    data(pp).Nch = sum(data(pp).rt_sel) + sum(data(pp).acc_sel);
    
end


%% === PARAMETERS =========================================================
np = 0;

clear paramName link reg



%  ============================  RT PARAMETERS  ===========================
%  rt_sigma  =====================
np = np + 1;
paramName{np}   = 'rt_{sigma}';
link{np}        = @(x) log(1+exp(x));
reg{np}         = [];



%  rt_ndt  =====================
np = np + 1;
paramName{np}   = 'rt_{ndt}';
link{np}        = @(x) (.5)./(1+exp(-x));
reg{np}         = [];


%  rt_1  =====================
np = np + 1;
paramName{np}   = 'rt_{1}';
link{np}        = @(x) x;
reg{np}         = [];


%  rt_prev  =====================
np = np + 1;
paramName{np}   = 'rt_{ar}';
link{np}        = @(x) x;
reg{np}         = [];



%  rt_targ  =====================
if ~strcmp(name, 'exp1')
    np = np + 1;
    paramName{np}   = 'rt_{targ}';
    link{np}        = @(x) x;
    reg{np}         = [];
end

%  rt_dist  =====================
np = np + 1;
paramName{np}   = 'rt_{dist}';
link{np}        = @(x) x;
reg{np}         = [];


%  rt_targDist  =====================
if ~strcmp(name, 'exp1')
    np = np + 1;
    paramName{np}   = 'rt_{targDist}';
    link{np}        = @(x) x;
    reg{np}         = [];
end





%  ============================  ACC PARAMETERS  ==========================
%  acc_s  =====================
np = np + 1;
paramName{np}   = 'acc_{lapse}';
link{np}        = @(x) (1)./(1+exp(-x));
reg{np}         = [];


%  acc_1  =====================
np = np + 1;
paramName{np}   = 'acc_{1}';
link{np}        = @(x) x;
reg{np}         = [];


%  acc_prev  =====================
np = np + 1;
paramName{np}   = 'acc_{ar}';
link{np}        = @(x) x;
reg{np}         = [];


%  acc_targ  =====================
np = np + 1;
paramName{np}   = 'acc_{targ}';
link{np}        = @(x) x;
reg{np}         = [];


%  acc_dist  =====================
np = np + 1;
paramName{np}   = 'acc_{dist}';
link{np}        = @(x) x;
reg{np}         = [];



%  acc_targ  =====================
if ~strcmp(name, 'exp1')
    np = np + 1;
    paramName{np}   = 'acc_{targDist}';
    link{np}        = @(x) x;
    reg{np}         = [];
end



%  pwr_{targ}  =====================
if ~strcmp(name, 'exp1')
    np = np + 1;
    paramName{np}   = 'pwr_{targ}';
    link{np}        = @(x) log(1+exp(x));
    reg{np}         = [];
end


%  pwr_{dist}  =====================
np = np + 1;
paramName{np}   = 'pwr_{dist}';
link{np}        = @(x) log(1+exp(x));
reg{np}         = [];








for tt = 1:length(data)
    
    data(tt).link = link;
    data(tt).paramName = paramName;
    data(tt).name = name;
    data(tt).likfun = likfun;
    
end




%% === FIT MODEL ==========================================================

if serial
    
    [E,V,alpha,stats,bf,fitparams] = semfit(likfun, data, np, reg,[],[],1,[],[],['./in-progress/', datestr(now, 'yyyy-mm-dd_HH-MM'), '_', saveAs]);
    
else
    
    tic;
    [E,V,alpha,stats,bf,fitparams] = emfit(likfun, data, np, reg,[],[],1,[],[],['./in-progress/', datestr(now, 'yyyy-mm-dd_HH-MM'), '_', saveAs]);
    
    fitDuration = toc/60
    
end


%% save results

% add fit info to results
results.model           = saveAs;
results.likfun          = likfun;

results.stats           = stats;
results.fitparams       = fitparams;

results.paramName       = paramName;
results.link            = link;

results.alpha           = alpha;
results.E               = E;
results.V               = V;
results.bf              = bf;

results.fitDuration     = fitDuration;


disp(' ')
disp(' parameters ')
disp(' ')
disp(paramName)
disp(' ')
disp(results.alpha')
disp(' ')



%% SAVE
save_fn = [datestr(datetime,'yyyy-mm-dd_HH-MM_'), saveAs]
save(['./fit-results/', save_fn], 'data', 'results');














end













