function [lik, out] = lik_rewMotionStatic(p, data, mu, nui, doprior, llopt)


%% ---------------------------------------------- get prior

lik = 0;

if doprior == 1 % if group prior
    
    [lik] = logGaussianPrior(p, mu, nui, doprior);
    
end

%% ---------------------------------------------- Load Parameters


% ===== transform parameters =====

for nn = 1:length(data.link)
    
    p(nn) = data.link{nn}(p(nn));
    
end



% ===== assign parameters =====

np = 1;



%  RT PARAMETERS  ===========================
rt_sigma    = p(np);np = np+1;

rt_ndt      = p(np);np = np+1;


% intercept
rt_1        = p(np);np = np+1;

% AR
rt_ar     = p(np);np = np+1;

% rew
rt_rew      = p(np);np = np+1;

% within-trial
rt_targ     = p(np);np = np+1;

rt_dist     = p(np);np = np+1;

% rew mod
rt_rewTarg  = p(np);np = np+1;

rt_rewDist     = p(np);np = np+1;






%  ACC PARAMETERS  ==========================
acc_lapse   = p(np);np = np+1;

% intercept
acc_1        = p(np);np = np+1;

% AR
acc_ar     = p(np);np = np+1;

% rew
acc_rew      = p(np);np = np+1;

% within-trial
acc_targ     = p(np);np = np+1;

acc_dist     = p(np);np = np+1;

% rew mod
acc_rewTarg  = p(np);np = np+1;

acc_rewDist  = p(np);np = np+1;



% power
if ~strcmp(data.name, 'exp1')
    pwr_targ = p(np);np = np+1;
else
    pwr_targ = 0;
end

pwr_dist = p(np);np = np+1;




%% ---------------------------------------------- Load Behavior
% center = @(x) x-nanmean(x);


rt = data.rt0;
prevRT = data.rt1;
rt_sel = data.rt_sel;


% acc = data.acc0;
acc = data.ptResp;
prevAcc = data.sgnPtResp1;
acc_sel = data.acc_sel;



% within-trial
if ~strcmp(data.name, 'exp1')
    targ = tanh(pwr_targ.*(data.targ0 + .6))./tanh(pwr_targ);
    respTarg = targ.*data.sgnResp;
else
    targ = 0*data.dist0;
    respTarg = data.sgnResp;
end


dist = tanh(pwr_dist*data.dist0)./tanh(pwr_dist);
respDist =  dist.*data.sgnResp;



% reward info
rew = 2*(grp2idx(data.rew) - nanmean(grp2idx(data.rew)));
rewTarg = rew.*targ;
rewDist = rew.*dist;

rewRespTarg = rew.*respTarg;
rewRespDist = rew.*respDist;



%% ---------------------------------------------- RT FIT

cTarg = targ(rt_sel) - mean(targ(rt_sel));
cPrevRt = prevRT(rt_sel) - mean(prevRT(rt_sel));

% organize
rt_beta     = [...
    rt_1, rt_ar, rt_targ, rt_dist, ...
    rt_rew, rt_rewTarg, rt_rewDist]';

rt_X        = [...
    ones(sum(rt_sel), 1), cPrevRt, cTarg, dist(rt_sel),...
    rew(rt_sel), cTarg.*rew(rt_sel), rewDist(rt_sel)];


% get pdf
rt_yhat     = rt_X * rt_beta;
rt_pdf      = lognpdf(rt(rt_sel)-rt_ndt, rt_yhat, rt_sigma);

rt_pdf(rt_pdf<realmin | ~isfinite(rt_pdf)) = realmin;


% get loglik
lik = lik + sum(log(rt_pdf));



%% ---------------------------------------------- ACC FIT


acc_lapse = 1 ./ (1+exp(-(acc_lapse + acc_rew*rew(acc_sel))));


% organize
acc_beta     = [...
    acc_1, acc_ar, acc_targ, acc_dist,...
    0, acc_rewTarg, acc_rewDist]';

acc_X        = [...
    ones(sum(acc_sel), 1), prevAcc(acc_sel), respTarg(acc_sel), respDist(acc_sel),...
    rew(acc_sel), rewRespTarg(acc_sel), rewRespDist(acc_sel)];


% get pdf
acc_yhat = (1-acc_lapse) ./ (1+exp(-acc_X*acc_beta)) + acc_lapse*.50;
acc_pdf = acc_yhat .* acc(acc_sel) + (1 - acc_yhat) .* (1 - acc(acc_sel));

acc_pdf(acc_pdf<realmin | ~isfinite(acc_pdf)) = realmin;


% get loglik
lik = lik + sum(log(acc_pdf));





%% =================================================

lik = -lik;




%% save variables
if nargout > 1
    
    % data
    out.data = data;
    
    
    % ===== RT
    out.rt.rt_beta = rt_beta;
    out.rt.rt_X = rt_X;
    out.rt.rt = rt;
    out.rt.rt_sel = rt_sel;
    
    out.rt.rt_yhat = rt_yhat;
    out.rt.rt_pdf = rt_pdf;
    
    
    % ===== ACC
    out.acc.acc_beta = acc_beta;
    out.acc.acc_X = acc_X;
    out.acc.acc = acc;
    out.acc.acc_sel = acc_sel;
    
    
    out.acc.acc_yhat = acc_yhat;
    out.acc.acc_pdf = acc_pdf;
    
    
end



