function [lik, out] = lik_withinStaticNoX(p, data, mu, nui, doprior, llopt)


%% ---------------------------------------------- get prior

lik = 0;

if doprior == 1 % if group prior
    
    [lik] = logGaussianPrior(p, mu, nui, doprior);
    
end


%% ---------------------------------------------- Fixed Parameters



%% ---------------------------------------------- Load Parameters


% ===== transform parameters =====

for nn = 1:length(data.link)
    try
        p(nn) = data.link{nn}(p(nn));
    catch
        keyboard
    end
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


if ~strcmp(data.name, 'exp1')
    rt_targDist     = 0;%p(np);np = np+1;
else
    rt_targDist = 0;
end





%  ACC PARAMETERS  ==========================
acc_lapse   = p(np);np = np+1;


% intercept
acc_1       = p(np);np = np+1;

% AR
acc_ar     = p(np);np = np+1;

% within-trial
acc_targ    = p(np);np = np+1;

acc_dist    = p(np);np = np+1;


if ~strcmp(data.name, 'exp1')
    acc_targDist    = 0;%p(np);np = np+1;
else
    acc_targDist = 0;
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





%% ---------------------------------------------- RT FIT

cTarg = targ(rt_sel) - mean(targ(rt_sel));
cPrevRt = prevRT(rt_sel) - mean(prevRT(rt_sel));

% organize
rt_beta     = [rt_1, rt_ar, rt_targ, rt_dist, rt_targDist]';
rt_X        = [ones(sum(rt_sel), 1), cPrevRt, cTarg, dist(rt_sel), cTarg.*dist(rt_sel)];


% get pdf
rt_yhat     = rt_X * rt_beta;
rt_pdf      = lognpdf(rt(rt_sel)-rt_ndt, rt_yhat, rt_sigma);

rt_pdf(rt_pdf<realmin | ~isfinite(rt_pdf)) = realmin;


% get loglik
lik = lik + sum(log(rt_pdf));




%% ---------------------------------------------- ACC FIT


% organize
acc_beta     = [acc_1, acc_ar, acc_targ, acc_dist, acc_targDist]';
acc_X        = [ones(sum(acc_sel), 1), prevAcc(acc_sel), respTarg(acc_sel), respDist(acc_sel), respTarg(acc_sel).*respDist(acc_sel)];


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



