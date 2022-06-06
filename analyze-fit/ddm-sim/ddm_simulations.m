%% dynamics sim for PACT



%% setup

% clear
clear; clc;

% add toolbox
addpath(genpath('./dm-util/'));

% load colors
load('batlow');
cols = batlow;
cIdx = round(linspace(1, length(cols), 7));
xoff = linspace(-.2,.2,5);

% functions
center = @(x) x - nanmean(x);

%% parameters


% sim parameters
nTrials     = 1e6;
dt          = 1e-3;
nDisc       = 5;
saveLoc     = ['./model_fits/', datestr(datetime,'yyyy-mm-dd_HH-MM')]

mkdir(saveLoc)
mkdir([saveLoc, '/fig'])
mkdir([saveLoc, '/png'])
mkdir([saveLoc, '/eps'])
mkdir([saveLoc, '/pdf'])


% task parameters
tr_dist = datasample(linspace(-1, 1, 7)', nTrials, 1);
tr_targ = datasample(linspace(.3, .9, 7)', nTrials, 1);


% accum param
targ_wt     = ones(nTrials,1) * 1.5;
dist_wt     = ones(nTrials,1) * .15;
bound       = ones(nTrials,1) * 1;
noise       = ones(nTrials,1) * 1;


% variability param
s_targ  = unifrnd(-.5, .5, nTrials, 1);
s_dist  = unifrnd(-.05, .05, nTrials, 1);
s_bound = unifrnd(-.25, .25, nTrials, 1);
s_noise = unifrnd(-.25, .25, nTrials, 1);



% models
vars = {'none',...
    'targ', 'dist', 'gain', 'bound', 'noise',...
    'none', 'none', 'none', 'none', 'none'}

dyns = {'none',...
    'none', 'none', 'none', 'none', 'none', ...
    'attn', 'attnBound', 'boundDown', 'noiseUp', 'noiseDown'}



assert(length(vars) == length(dyns))




%% run simulations

clear mdl ffx
for ss = 1:length(vars)

    % between-trial variability

    switch vars{ss}

        case 'none'

            mu = dist_wt.*tr_dist + targ_wt.*tr_targ;
            bound = bound;
            noise = noise;

        case 'targ'

            mu = dist_wt.*tr_dist + (targ_wt + s_targ).*tr_targ;
            bound = bound;
            noise = noise;

        case 'dist'

            mu = (dist_wt + s_dist).*tr_dist + targ_wt.*tr_targ;
            bound = bound;
            noise = noise;

        case 'gain'

            mu = (dist_wt + s_dist).*tr_dist + (targ_wt + s_targ).*tr_targ;
            bound = bound;
            noise = noise;

        case 'bound'

            mu = dist_wt.*tr_dist + targ_wt.*tr_targ;
            bound = bound + s_bound;
            noise = noise;

        case 'noise'

            mu = (dist_wt.*tr_dist + targ_wt.*tr_targ);
            bound = bound;
            noise = noise + s_noise;


    end



    fprintf('\n===== var: %s | dyn: %s ===== \n\n', vars{ss}, dyns{ss})


    % sim
    tic
    [tr_rt, tr_acc] = deal(nan(nTrials,1));
    parfor tt = 1:nTrials


        switch dyns{ss}

            case 'none'

                [tr_rt(tt), tr_acc(tt)] = ddm_rand_full(...
                    mu(tt),...
                    noise(tt),...
                    -bound(tt),...
                    bound(tt),...
                    dt, 1);


            case 'attn'

                [tr_rt(tt), tr_acc(tt)] = ddm_rand_full(...
                    linspace(.2, 0, 1/dt).*tr_dist(tt) + linspace(1, 3, 1/dt).*tr_targ(tt),...
                    noise(tt),...
                    -bound(tt),...
                    bound(tt),...
                    dt, 1);


            case 'attnBound'

                [tr_rt(tt), tr_acc(tt)] = ddm_rand_full(...
                    linspace(.2, 0, 1/dt).*tr_dist(tt) + linspace(1, 3, 1/dt).*tr_targ(tt),...
                    noise(tt),...
                    -linspace(1.5, .1, 2/dt),...
                    linspace(1.5, .1, 2/dt),...
                    dt, 1);


            case 'boundDown'


                [tr_rt(tt), tr_acc(tt)] = ddm_rand_full(...
                    mu(tt),...
                    noise(tt),...
                    -linspace(1.5, .1, 2/dt),...
                    linspace(1.5, .1, 2/dt),...
                    dt, 1);


            case 'noiseUp'


                [tr_rt(tt), tr_acc(tt)] = ddm_rand_full(...
                    mu(tt),...
                    linspace(1, 2, 1/dt),...
                    -bound(tt),...
                    bound(tt),...
                    dt, 1);


            case 'noiseDown'


                [tr_rt(tt), tr_acc(tt)] = ddm_rand_full(...
                    mu(tt),...
                    linspace(1, .5, 1/dt),...
                    -bound(tt),...
                    bound(tt),...
                    dt, 1);


        end





    end
    toc

    meanRT  = mean(tr_rt)
    stdRT  = std(tr_rt)
    meanACC = mean(tr_acc)

    % fit data w/ regression
    tbl = table;
    tbl.acc = tr_acc;
    tbl.crt = tr_rt - nanmean(tr_rt);
    tbl.targ = tr_targ;
    tbl.dist = tr_dist;

    % regression
    mdl = fitglm(tbl, 'acc ~ (dist + targ)*crt', 'Distribution', 'binomial');
    ffx{ss} = mdl.Coefficients;



    % plot dynamics ====================
    tbl.discRT = discretize(tiedrank(tbl.crt), nDisc);
    tbl.discTarg = discretize(tbl.targ, 7);
    tbl.discDist = discretize(tbl.dist, 7);


    % start figure
    f=figure('units','inch','position',[0,0,7,5]); hold on;
    tiledlayout(1,2, 'TileSpacing','compact', 'Padding','compact');



    % === target
    gt = grpstats(tbl, {'discRT', 'discTarg'});

    nexttile; hold on;
    title(sprintf('targ || var: %s, dyn: %s', vars{ss}, dyns{ss}))
    set(gca, 'LineWidth', 1, 'TickDir', 'out')


    for ds = 1:nDisc

        plot(gt.discTarg(gt.discRT == ds) + xoff(ds), gt.mean_acc(gt.discRT == ds),...
            '-o', 'LineWidth', 2,...
            'color', cols(cIdx(ds+1),:), 'MarkerFaceColor', cols(cIdx(ds+1),:))

    end
    xticks([])
    xlabel('Target Coherence')
    ylabel('accuarcy')




    % === distractor
    gt = grpstats(tbl, {'discRT', 'discDist'});

    nexttile; hold on;
    title(sprintf('dist || var: %s, dyn: %s', vars{ss}, dyns{ss}))
    set(gca, 'LineWidth', 1, 'TickDir', 'out')
    yline(0, 'LineWidth', 1)


    for ds = 1:nDisc

        plot(gt.discDist(gt.discRT == ds) + xoff(ds), center(gt.mean_acc(gt.discRT == ds)),...
            '-o', 'LineWidth', 2, ...
            'color', cols(cIdx(ds+1),:), 'MarkerFaceColor', cols(cIdx(ds+1),:))

    end
    xticks([])
    xlabel('Distractor Congruence')
    ylabel('accuarcy (RT centered)')


    % save figures
    saveas(f, sprintf('%s/fig/var-%s_dyn-%s', saveLoc, vars{ss}, dyns{ss}), 'fig')
    saveas(f, sprintf('%s/png/var-%s_dyn-%s', saveLoc, vars{ss}, dyns{ss}), 'png')
    saveas(f, sprintf('%s/pdf/var-%s_dyn-%s', saveLoc, vars{ss}, dyns{ss}), 'pdf')
    close(f)

    % save model
    save([saveLoc, '/fit.mat'], 'ffx', 'vars', 'dyns', '-v7.3')


end



fprintf('\ndone!\n\n')




%% display fit
clc;
disp('standard'); ffx{ismember(vars, 'none') & ismember(dyns, 'none')}
disp('attn'); ffx{ismember(vars, 'none') & ismember(dyns, 'attn')}
disp('attnbound'); ffx{ismember(vars, 'none') & ismember(dyns, 'attnBound')}

disp('targVar'); ffx{ismember(vars, 'targ') & ismember(dyns, 'none')}
disp('distVar'); ffx{ismember(vars, 'dist') & ismember(dyns, 'none')}
disp('noise'); ffx{ismember(vars, 'noise') & ismember(dyns, 'none')}
disp('bound'); ffx{ismember(vars, 'bound') & ismember(dyns, 'none')}

disp('noiseUp'); ffx{ismember(vars, 'none') & ismember(dyns, 'noiseUp')}
disp('noiseDown'); ffx{ismember(vars, 'none') & ismember(dyns, 'noiseDown')}
disp('boundDown'); ffx{ismember(vars, 'none') & ismember(dyns, 'boundDown')}





