function [lik, out] = lik_adaptStatic(p, data, mu, nui, doprior, llopt)


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


% within-trial
if ~strcmp(data.name, 'exp1')
    rt_targ     = p(np);np = np+1;
else
    rt_targ = 0;
end

rt_dist     = p(np);np = np+1;


% previous-trial
if ~strcmp(data.name, 'exp1')
    rt_targ1     = p(np);np = np+1;
else
    rt_targ1 = 0;
end

rt_dist1     = p(np);np = np+1;


% adaptation
rt_distDist1     = p(np);np = np+1;

if ~strcmp(data.name, 'exp1')
    rt_distTarg1     = p(np);np = np+1;
    rt_targDist1     = p(np);np = np+1;
    rt_targTarg1     = p(np);np = np+1;
else
    rt_distTarg1     = 0;
    rt_targDist1     = 0;
    rt_targTarg1     = 0;
end





%  ACC PARAMETERS  ==========================
acc_lapse   = p(np);np = np+1;


% intercept
acc_1       = p(np);np = np+1;

% AR
acc_ar      = p(np);np = np+1;

% within-trial
acc_targ     = p(np);np = np+1;


acc_dist     = p(np);np = np+1;


% previous-trial
if ~strcmp(data.name, 'exp1')
    acc_targ1     = p(np);np = np+1;
else
    acc_targ1 = 0;
end

acc_dist1     = p(np);np = np+1;




% adaptation
acc_distDist1     = p(np);np = np+1;

if ~strcmp(data.name, 'exp1')
    acc_distTarg1     = p(np);np = np+1;
    acc_targDist1     = p(np);np = np+1;
    acc_targTarg1     = p(np);np = np+1;
else
    acc_distTarg1     = 0;
    acc_targDist1     = 0;
    acc_targTarg1     = 0;
end



% power
if ~strcmp(data.name, 'exp1')
    pwr_targ = p(np);np = np+1;
else
    pwr_targ = 0;
end

pwr_dist = p(np);np = np+1;




%% ---------------------------------------------- Load Behavior

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



% previous trial
if ~strcmp(data.name, 'exp1')
    %     targ1 = log(1 + (data.targ1 + .6));
    targ1 = data.targ1;
else
    targ1 = 0*data.dist1;
end

dist1 = data.dist1;



% interactions

distDist1 = dist.*dist1;
respDistDist1 = respDist.*dist1;

if ~strcmp(data.name, 'exp1')
    
    distTarg1     = dist.*targ1;
    targDist1     = targ.*dist1;
    targTarg1     = targ.*targ1;
    
else
    
    distTarg1     = dist*0;
    targDist1     = dist*0;
    targTarg1     = dist*0;
    
    respDistTarg1     = dist*0;
    respTargDist1     = dist*0;
    respTargTarg1     = dist*0;
    
end





%% ---------------------------------------------- RT FIT

cTarg = targ(rt_sel) - mean(targ(rt_sel));
cTarg1 = targ1(rt_sel) - mean(targ1(rt_sel));
cPrevRt = prevRT(rt_sel) - mean(prevRT(rt_sel));

% organize
rt_beta     = [...
    rt_1, rt_ar, rt_targ, rt_dist, rt_targ1, rt_dist1,...
    rt_distDist1, rt_distTarg1, rt_targDist1, rt_targTarg1]';

rt_X        = [...
    ones(sum(rt_sel), 1), cPrevRt, cTarg, dist(rt_sel), cTarg1, dist1(rt_sel),...
    distDist1(rt_sel), distTarg1(rt_sel), cTarg.*dist1(rt_sel), cTarg.*cTarg1];


% get pdf
rt_yhat     = rt_X * rt_beta;
rt_pdf      = lognpdf(rt(rt_sel)-rt_ndt, rt_yhat, rt_sigma);

rt_pdf(rt_pdf<realmin | ~isfinite(rt_pdf)) = realmin;


% get loglik
lik = lik + sum(log(rt_pdf));



%% ---------------------------------------------- ACC FIT
cTarg1 = targ1 - mean(targ1(acc_sel));

respDistTarg1     = respDist.*cTarg1;
respTargDist1     = respTarg.*dist1;
respTargTarg1     = respTarg.*cTarg1;

acc_lapse = 1 ./ (1+exp(-(acc_lapse + acc_targ1*targ1(acc_sel) + acc_dist1*dist1(acc_sel))));

% organize
acc_beta     = [...
    acc_1, acc_ar, acc_targ, acc_dist, 0, 0,...
    acc_distDist1, acc_distTarg1, acc_targDist1, acc_targTarg1]';

acc_X        = [...
    ones(sum(acc_sel), 1), prevAcc(acc_sel), respTarg(acc_sel), respDist(acc_sel), cTarg1(acc_sel), dist1(acc_sel),...
    respDistDist1(acc_sel), respDistTarg1(acc_sel), respTargDist1(acc_sel), respTargTarg1(acc_sel)];


% get pdf
acc_yhat = (1 - acc_lapse) ./ (1+exp(-acc_X*acc_beta)) + acc_lapse*0.5;
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



