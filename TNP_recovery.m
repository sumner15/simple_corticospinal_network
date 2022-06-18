%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% purpose: A firing-rate neural network model of M1 and SMA before and
% after stroke. In addition, this model is further driven to simulate
% traditional therapy as well as targeted feedback.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc; set(0,'defaultlinelinewidth',2.5)

%% set parameters here
N = 1000;               %number of neurons
maxRate = 100;          %max firing rate
maxStroke = round(N/3); %max stroke size (in # of affected cells)
target = [];            %updated under targeted feedback simulation
chem = false;           %whether or not to use chemotaxis
nDays = 180;            %number of days to run the simulation
dose_days = 20000;      % number of trials in f(trials) sim
dose_trials = ones(1,dose_days); 
video = false;          %true=create and save videos
nMC = 5;                %number of monte carlo sims for shaded error

%% network setup
% multiplier for weak/strong/sluggish/variable neurons
mult = 1.5;
% neuron weightings: WEIGHTING (w) = log normal dist. of (mean,var,N)
w = lognrndWrap(1,.25,N);

weakI = 1:round(N/10);          %weakly connected 
strongI = (weakI(end)+1):N;     %strongly connected
w(weakI) = w(weakI)/mult;          
w(strongI) = w(strongI)*mult;   
weightGroup(weakI) = {'weak connection'};
weightGroup(strongI) = {'strong connection'};

%stochastic standard deviation (fast/slow circuits)
ssd = lognrndWrap(1,.2,N);

[~,varI] = datasample(ssd,N/2,'Replace',false);     %variable indices
sluI = 1:N; sluI(varI) = [];                        %sluggish indices
ssd(varI) = ssd(varI)*mult;
ssd(sluI) = ssd(sluI)/mult;
speedGroup(1:N) = {'sluggish neuron'};
speedGroup(varI) = {'variable neuron'};

%save for later (they may be altered)
w0 = w; ssd0 = ssd;             

%% initialize simulation
% log normal distribution of initial firing rates 
v = maxRate/8; % variance
m = maxRate/4; % mean
X0 = lognrndWrap(m,v,N);
X00 = X0;      % save for later

% running simulation as f(trials)
for i = 1:nMC
    [XfInit_trial, FInit_trial(i,:), fMaxInit_trial] = ...
        simulate_CS(X0,w,maxRate,ssd,target,chem,dose_days,2*dose_trials);
end

%% injury simulation

% running simulation as f(trials)
for i = 1:nMC
    % make random stroke
    strokeInds = randperm(N, maxStroke);
    w_stroke = w;
    w_stroke(strokeInds) = 0;
    X0 = X00;
    X0(strokeInds) = 0;
    % run sim
    [XfStroke_trial, FStroke_trial(i,:), fMaxStroke_trial] = ...
        simulate_CS(X0,w_stroke,maxRate,ssd,target,chem,dose_days,dose_trials);
end

%% targeted feedback simulation

% target = sort(intersect(strongI,sluI)); targetStr = 'sluggish()strong';
target = sort(sluI); targetStr = 'all sluggish';

% running simulation as f(trials)
for i = 1:nMC
    % make random stroke
    strokeInds = randperm(N, maxStroke);
    w_stroke = w;
    w_stroke(strokeInds) = 0;
    X0 = X00;
    X0(strokeInds) = 0;
    % run sim
    [XfTarget_trial, FTarget_trial(i,:), fMaxTarget_trial] = ...
        simulate_CS(X0,w_stroke,maxRate,ssd,target,chem,dose_days,dose_trials);
end

%% save results for later use and begin plotting
save(strcat('./results/TNP_recovery_result_',datestr(now,'yyyy-mm-dd')));
plotResults()
