function [lik, out] = lik_adaptColMotStatic(p, data, mu, nui, doprior, llopt)


%% ---------------------------------------------- get prior

lik = 0;

if doprior == 1 % if group prior
    
    [lik] = logGaussianPrior(p, mu, nui, doprior);
    
end



%% ---------------------------------------------- Fixed Parameters



%% ---------------------------------------------- Load Parameters


% ===== transform parameters =====

% for nn = 1:length(data.link)
%     try
%         p(nn) = data.link{nn}(p(nn));
%     catch
%         keyboard
%     end
% end



% ===== assign parameters =====

np = 1;



% ================= COLOR BLOCKS


%  RT PARAMETERS  ===========================
rt_sigma    = data.link{np}(p(np)); np = np+1;

rt_ndt      = data.link{np}(p(np)); np = np+1;

% intercept
rt_1        = data.link{np}(p(np)); np = np+1;

rt_ar        = data.link{np}(p(np)); np = np+1;


% within-trial
if ~strcmp(data.name, 'exp1')
    rt_targ     = data.link{np}(p(np)); np = np+1;
else
    rt_targ = 0;
end

rt_dist     = data.link{np}(p(np)); np = np+1;


% adapt simple
if ~strcmp(data.name, 'exp1')
    rt_targ1     = data.link{np}(p(np)); np = np+1;
else
    rt_targ1 = 0;
end
rt_dist1        = data.link{np}(p(np)); np = np+1;


% adapt x
if ~strcmp(data.name, 'exp1')
    rt_targTarg1     = data.link{np}(p(np)); np = np+1;
else
    rt_targTarg1 = 0;
end

npdd_rt = np;
rt_distDist1     = data.link{np}(p(np)); np = np+1;




%  ACC PARAMETERS  ==========================
acc_lapse   = data.link{np}(p(np)); np = np+1;


% intercept
acc_1       = data.link{np}(p(np)); np = np+1;

acc_ar      = data.link{np}(p(np)); np = np+1;

% within-trial
acc_targ    = data.link{np}(p(np)); np = np+1;

acc_dist    = data.link{np}(p(np)); np = np+1;


% adapt simple
if ~strcmp(data.name, 'exp1')
    acc_targ1     = data.link{np}(p(np)); np = np+1;
else
    acc_targ1 = 0;
end

acc_dist1    = data.link{np}(p(np)); np = np+1;


% adapt x
if ~strcmp(data.name, 'exp1')
    acc_targTarg1     = data.link{np}(p(np)); np = np+1;
else
    acc_targTarg1 = 0;
end

npdd_acc = np;
acc_distDist1    = data.link{np}(p(np)); np = np+1;






% power
if ~strcmp(data.name, 'exp1')
    pwr_targ = data.link{np}(p(np)); np = np+1;
else
    pwr_targ = 0;
end

pwr_dist = data.link{np}(p(np)); np = np+1;





% ================= MOTION BLOCKS

mp = 1;

%  RT PARAMETERS  ===========================
rtMot_sigma    = data.link{np}(p(mp) + p(np)); np = np+1; mp = mp+1;

rtMot_ndt      = data.link{np}(p(mp) + p(np)); np = np+1; mp = mp+1;

% intercept
rtMot_1        = data.link{np}(p(mp) + p(np)); np = np+1; mp = mp+1;

rtMot_ar        = data.link{np}(p(mp) + p(np)); np = np+1; mp = mp+1;


% within-trial
if ~strcmp(data.name, 'exp1')
    rtMot_targ     = data.link{np}(p(mp) + p(np)); np = np+1; mp = mp+1;
else
    rtMot_targ = 0;
end

rtMot_dist     = data.link{np}(p(mp) + p(np)); np = np+1; mp = mp+1;


% adapt simple
if ~strcmp(data.name, 'exp1')
    rtMot_targ1     = data.link{np}(p(mp) + p(np)); np = np+1; mp = mp+1;
else
    rtMot_targ1 = 0;
end

rtMot_dist1     = data.link{np}(p(mp) + p(np)); np = np+1; mp = mp+1;


% adapt x
if ~strcmp(data.name, 'exp1')
    rtMot_targTarg1     = data.link{np}(p(mp) + p(npdd_rt)); np = np+1; mp = mp+1;
else
    rtMot_targTarg1 = 0;
end

rtMot_distDist1     = data.link{np}(p(mp) + p(npdd_rt)); np = np+1; mp = mp+1;



%  ACC PARAMETERS  ==========================
accMot_lapse   = data.link{np}(p(mp) + p(np)); np = np+1; mp = mp+1;

% intercept
accMot_1       = data.link{np}(p(mp) + p(np)); np = np+1; mp = mp+1;

accMot_ar       = data.link{np}(p(mp) + p(np)); np = np+1; mp = mp+1;

% within-trial
accMot_targ    = data.link{np}(p(mp) + p(np)); np = np+1; mp = mp+1;

accMot_dist    = data.link{np}(p(mp) + p(np)); np = np+1; mp = mp+1;


% adapt simple
if ~strcmp(data.name, 'exp1')
    accMot_targ1     = data.link{np}(p(mp) + p(np)); np = np+1; mp = mp+1;
else
    accMot_targ1 = 0;
end

accMot_dist1    = data.link{np}(p(mp) + p(np)); np = np+1; mp = mp+1;


% adapt x
if ~strcmp(data.name, 'exp1')
    accMot_targTarg1     = data.link{np}(p(mp) + p(npdd_acc)); np = np+1; mp = mp+1;
else
    accMot_targTarg1 = 0;
end

accMot_distDist1    = data.link{np}(p(mp) + p(npdd_acc)); np = np+1; mp = mp+1;



%% ---------------------------------------------- Load Behavior

rt = data.rt0;
prevRT = data.rt1;
rt_sel = data.rt_sel;
rtMot_sel = data.rtMot_sel;
% rt_sel = data.rtMot_sel;

% acc = data.acc0;
acc = data.ptResp;
prevAcc = data.sgnPtResp1;
acc_sel = data.acc_sel;
accMot_sel = data.accMot_sel;
% acc_sel = data.accMot_sel;

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
    rt_distDist1, 0, 0, rt_targTarg1]';

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
    acc_distDist1, 0, 0, acc_targTarg1]';

acc_X        = [...
    ones(sum(acc_sel), 1), prevAcc(acc_sel), respTarg(acc_sel), respDist(acc_sel), cTarg1(acc_sel), dist1(acc_sel),...
    respDistDist1(acc_sel), respDistTarg1(acc_sel), respTargDist1(acc_sel), respTargTarg1(acc_sel)];


% get pdf
acc_yhat = (1 - acc_lapse) ./ (1+exp(-acc_X*acc_beta)) + acc_lapse*0.5;
acc_pdf = acc_yhat .* acc(acc_sel) + (1 - acc_yhat) .* (1 - acc(acc_sel));

acc_pdf(acc_pdf<realmin | ~isfinite(acc_pdf)) = realmin;



% get loglik
lik = lik + sum(log(acc_pdf));










%% ---------------------------------------------- RT FIT

cTarg = targ(rtMot_sel) - mean(targ(rtMot_sel));
cTarg1 = targ1(rtMot_sel) - mean(targ1(rtMot_sel));
cPrevRt = prevRT(rtMot_sel) - mean(prevRT(rtMot_sel));

% organize
rtMot_beta     = [...
    rtMot_1, rtMot_ar, rtMot_targ, rtMot_dist, rtMot_targ1, rtMot_dist1,...
    rtMot_distDist1, 0, 0, rtMot_targTarg1]';

rtMot_X        = [...
    ones(sum(rtMot_sel), 1), cPrevRt, cTarg, dist(rtMot_sel), cTarg1, dist1(rtMot_sel),...
    distDist1(rtMot_sel), distTarg1(rtMot_sel), cTarg.*dist1(rtMot_sel), cTarg.*cTarg1];


% get pdf
rtMot_yhat     = rtMot_X * rtMot_beta;
rtMot_pdf      = lognpdf(rt(rtMot_sel)-rtMot_ndt, rtMot_yhat, rtMot_sigma);

rtMot_pdf(rtMot_pdf<realmin | ~isfinite(rtMot_pdf)) = realmin;


% get loglik
lik = lik + sum(log(rtMot_pdf));



%% ---------------------------------------------- ACC FIT
cTarg1 = targ1 - mean(targ1(accMot_sel));

respDistTarg1     = respDist.*cTarg1;
respTargDist1     = respTarg.*dist1;
respTargTarg1     = respTarg.*cTarg1;

accMot_lapse = 1 ./ (1+exp(-(accMot_lapse + acc_targ1*targ1(accMot_sel) + acc_dist1*dist1(accMot_sel))));

% organize
accMot_beta     = [...
    accMot_1, accMot_ar, accMot_targ, accMot_dist, 0, 0,...
    accMot_distDist1, 0, 0, accMot_targTarg1]';

accMot_X        = [...
    ones(sum(accMot_sel), 1), prevAcc(accMot_sel), respTarg(accMot_sel), respDist(accMot_sel), cTarg1(accMot_sel), dist1(accMot_sel),...
    respDistDist1(accMot_sel), respDistTarg1(accMot_sel), respTargDist1(accMot_sel), respTargTarg1(accMot_sel)];


% get pdf
accMot_yhat = (1 - accMot_lapse) ./ (1+exp(-accMot_X*accMot_beta)) + accMot_lapse*0.5;
accMot_pdf = accMot_yhat .* acc(accMot_sel) + (1 - accMot_yhat) .* (1 - acc(accMot_sel));

accMot_pdf(accMot_pdf<realmin | ~isfinite(accMot_pdf)) = realmin;



% get loglik
lik = lik + sum(log(accMot_pdf));





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
    %     out.rt.rtMot_sel = rtMot_sel;
    
    out.rt.rt_yhat = rt_yhat;
    out.rt.rt_pdf = rt_pdf;
    
    
    % ===== ACC
    out.acc.acc_beta = acc_beta;
    out.acc.acc_X = acc_X;
    out.acc.acc = acc;
    out.acc.acc_sel = acc_sel;
    %     out.acc.accMot_sel = accMot_sel;
    
    
    out.acc.acc_yhat = acc_yhat;
    out.acc.acc_pdf = acc_pdf;
    
end



