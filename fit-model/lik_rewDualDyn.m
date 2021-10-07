function [lik, out] = lik_rewDualDyn(p, data, mu, nui, doprior, llopt)


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
rt_1       = p(np);np = np+1;

% previous
rt_ar       = p(np);np = np+1;

% rew
rt_rew       = p(np);np = np+1;

% accuracy
rt_acc    = p(np);np = np+1;

% rew
rt_rewAcc    = p(np);np = np+1;

% within-trial
rt_targ    = p(np);np = np+1;
rt_dist    = p(np);np = np+1;

% rew mod
rt_rewTarg    = p(np);np = np+1;
rt_rewDist    = p(np);np = np+1;

% acc mod
rt_targAcc = p(np);np = np+1;
rt_distAcc = p(np);np = np+1;

% rew x acc mod
rt_rewTargAcc = p(np);np = np+1;
rt_rewDistAcc = p(np);np = np+1;


%  ACC PARAMETERS  ==========================
acc_lapse   = p(np);np = np+1;

% intercept
acc_1        = p(np);np = np+1;

% AR
acc_ar      = p(np);np = np+1;

% rew
acc_rew      = p(np);np = np+1;

% crt
acc_crt      = p(np);np = np+1;

% rew crt
acc_rewCrt   = p(np);np = np+1;

% within-trial
acc_targ     = p(np);np = np+1;
acc_dist     = p(np);np = np+1;

% rew mod
acc_rewTarg     = p(np);np = np+1;
acc_rewDist     = p(np);np = np+1;

% crt mod
acc_crtTarg     = p(np);np = np+1;
acc_crtDist     = p(np);np = np+1;

% rew x crt mod
acc_rewCrtTarg     = p(np);np = np+1;
acc_rewCrtDist     = p(np);np = np+1;


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

correct = data.acc0;
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



% rew info
rew = 2*(grp2idx(data.rew) - nanmean(grp2idx(data.rew)));





%% ---------------------------------------------- RT FIT
crew    = rew(rt_sel);
cTarg = targ(rt_sel) - mean(targ(rt_sel));
cPrevRt = prevRT(rt_sel) - mean(prevRT(rt_sel));
cCorrect = correct(rt_sel) - mean(correct(rt_sel));

% organize
rt_beta     = [...
    rt_1, rt_ar, rt_targ, rt_dist, ...
    rt_rew, rt_acc, rt_rewAcc,...
    rt_rewTarg, rt_rewDist,...
    rt_targAcc, rt_distAcc,...
    rt_rewTargAcc, rt_rewDistAcc,...
    ]';

rt_X        = [...
    ones(sum(rt_sel), 1), cPrevRt, cTarg, dist(rt_sel),...
    crew, cCorrect, crew.*cCorrect,...
    cTarg.*crew, dist(rt_sel).*crew,...
    cTarg.*cCorrect, dist(rt_sel).*cCorrect,...
    cTarg.*cCorrect.*crew, dist(rt_sel).*cCorrect.*crew,...
    ];


% get pdf
rt_yhat     = rt_X * rt_beta;
rt_pdf      = lognpdf(rt(rt_sel)-rt_ndt, rt_yhat, rt_sigma);

rt_pdf(rt_pdf<realmin | ~isfinite(rt_pdf)) = realmin;


% get loglik
lik = lik + sum(log(rt_pdf));





%% ---------------------------------------------- ACC FIT
crew    = rew(acc_sel);
crt = rt(acc_sel) - nanmean(rt(acc_sel));
crt = crt - (crew\crt)*crew; %regress out reward effects from crt



acc_lapse = 1 ./ (1+exp(-(acc_lapse + acc_rew*crew + acc_crt*crt + acc_rewCrt.*crew.*crt)));


% organize
acc_beta     = [...
    acc_1, acc_ar,...
    acc_targ, acc_dist,...
    acc_rewTarg, acc_rewDist,...
    acc_crtTarg, acc_crtDist,...
    acc_rewCrtTarg, acc_rewCrtDist]';

acc_X        = [...
    ones(sum(acc_sel), 1), prevAcc(acc_sel),...
    respTarg(acc_sel), respDist(acc_sel),...
    crew.*respTarg(acc_sel), crew.*respDist(acc_sel),...
    crt.*respTarg(acc_sel), crt.*respDist(acc_sel),...
    crew.*crt.*respTarg(acc_sel), crew.*crt.*respDist(acc_sel),...
    ];


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
    out.acc.crt = crt;
    
    
    out.acc.acc_yhat = acc_yhat;
    out.acc.acc_pdf = acc_pdf;
    
    
end



