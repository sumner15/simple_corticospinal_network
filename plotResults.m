%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author: Sumner Norman
% script: plotResults
%
% Plots the network results for TNP_recovery
% takes no inputs, leaves no prisoners... I mean outputs
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% settings / which plots to produce
structurePlot = true;  % neuron parameters visualized
forcePlot = true;       % force production plots
covariancePlot = false; % covariance of force/activation sorted 

close all

%% load data and set common vars
while ~exist('FInit_trial','var')
   [filename, pathname]= uigetfile('*.mat','Select modeldata .mat file','./results');
   load(fullfile(pathname, filename));
end

legends = {'Before Learning', 'Scenario 1: Undamaged',...
    'Scenario 2: Conventional Rehab','Scenario 3: Targeted Rehab'};

brainInds{1} = intersect(strongI,varI);     %M1 contralateral
brainInds{2} = intersect(weakI,varI);       %M1 ipsilateral
brainInds{3} = intersect(strongI,sluI);     %SMA contralateral
brainInds{4} = intersect(weakI,sluI);       %SMA ipsilateral

colors{1} = [255 210 0]/285;    % dark yellow   (M1 cont)
colors{2} = [0 100 164]/285;    % dark blue     (M1 ips)
colors{3} = [224 225 95]/225;   % bright yellow (SMA cont)
colors{4} = [106 162 184]/200;  % bright blue   (SMA ips)

%% scatter/hist of neuron parameters
if structurePlot
    set(figure,'Position',[100 540 900 350])
    textY = 3150;
          

    subplot(121);       % HEMISPHERE ------------------------------
    histogram(w0(strongI),30,'FaceColor',colors{1})
    hold on
    histogram(w0(weakI),30,'FaceColor','blue')
    axis([0.15 10 0 N/6]); xlabel('Synaptic Weighting'); ylabel('# Neurons')
    h1_legend = legend('contralateral','ipsilateral');
    
    set(gca,'xscale','log'); grid on
    New_XTickLabel = get(gca,'xtick');
    set(gca,'XTickLabel',New_XTickLabel);
    
    text(0.5,textY,'Hemisphere','FontSize',16,'FontName','Times')
    
    subplot(122);       % MOTOR AREAS ------------------------------
    histogram(ssd0(varI),30,'FaceColor','black','FaceAlpha',0.8)
    hold on
    histogram(ssd0(sluI),30,'FaceColor','white')
    axis([0.15 10 0 N/6]); xlabel('Neuron Variability');
    h2_legend = legend('Primary','Secondary');
    
    set(gca,'xscale','log'); grid on
    New_XTickLabel = get(gca,'xtick');
    set(gca,'XTickLabel',New_XTickLabel);
    
    text(0.5,textY,'Motor Areas','FontSize',16,'FontName','Times')
end

%% plot simulation force results as f(trial)
if forcePlot     
    % subsample to make figures lighter
    nTrials = length(FInit_trial);
    t_plot = [1:100:nTrials nTrials];        
    % set useful variables   
    toPerc = 100/max(FInit_trial(:));  
    alphaVal = 0.4;   
        
    % plot force curves, scenario 1
    set(figure,'Position',[300 540 500 400])        
    [m,c] = getCI(FInit_trial*toPerc);
    shadedErrorBar(-1*fliplr(t_plot), m(t_plot), c(t_plot));
    hold on        
    % plot scenario 2
    [m,c] = getCI(FStroke_trial*toPerc);
    shadedErrorBar(t_plot, m(t_plot), c(t_plot), 'lineprops',{'markerfacecolor',colors{3}});
    % plot scenario 3
    [m,c] = getCI(FTarget_trial*toPerc);
    shadedErrorBar(t_plot, m(t_plot), c(t_plot), 'lineprops',{'markerfacecolor',colors{2}});
    % plot vertical line "curve" beetween/before after injury
    i1 = plot([0 1],[mean(FInit_trial(:,end)) FStroke_trial(1,1)]*toPerc,'-k');%injury line
    i1.Color(4) = alphaVal*2;   
    % plot maximum lines
    plot([-t_plot(end) -t_plot(1)],[100 100],':k')      %max before stroke
    plot([t_plot(1) t_plot(end)],[fMaxStroke_trial(end)*toPerc fMaxStroke_trial(end)*toPerc],':k') %injury line    
    text(nTrials/5,fMaxStroke_trial(1)*toPerc+3,'max recovery','fontsize',16)
    
    % set axes
    axis([-nTrials nTrials 20 102])    
    grid on   
    % set type
    xlabel('movement attempt'); ylabel('force (% of maximum)')    
    set(findall(gcf,'-property','FontSize'),'FontSize',20)          
    xticks(-nTrials:nTrials/2:nTrials)
    xticklabels({'0','10,000','20,000','30,000','40,000'})                
    % set legend    
    leg = legend(legends(2:4),'Location','SOUTHWEST');
    set(leg,'FontSize',12);
end

% stats
dist_target = FTarget_trial(:,end);
dist_stroke = FStroke_trial(:,end);
latent_recovery = (mean(dist_target) - mean(dist_stroke))...
    /(fMaxStroke_trial(end)-mean(dist_stroke));
[H, P, CI, stats] = ttest(dist_target, dist_stroke);


% printing results to console
fprintf('TNP reached %2.1f%% of the uninjured max\n', mean(FTarget_trial(:,end)*toPerc))
fprintf('That''s an increase of %2.1f%% compared to Scenario 2 (%2.1f%%).\n',...
    (mean(FTarget_trial(:,end))-mean(FStroke_trial(:,end)))*toPerc, mean(FStroke_trial(:,end)*toPerc))
fprintf('%2.1f%% recovery of latent capacity!\n',latent_recovery*100);
fprintf('t-test: t = %2.2f, p=%f\n', stats.tstat, P)


%% plot of covariance of force/activation sorted by cell parameters
if covariancePlot
    % calculate covariance/correlation coefficient/percent change
    covFX = corrcoef([FTarget' XHistory]); 
    warning('Is this what you expected? Could be wrong F...')
    covFX = covFX(2:end,1);
    percInc = 100*(XHistory(180,:)-XHistory(1,:))./XHistory(1,:);

    % sort weighting and cell variability
    [wSort, wSortI] = sort(w0);
    [ssdSort, ssdSortI] = sort(ssd0);

    % plotting
    set(figure,'Position',[75 25 1400 800])
    subplot(2,5,1:2)
    plot(wSort,covFX(wSortI),'.')
    xlim([.1 10])
    set(gca,'xscale','log'); grid on
    xlabel('weighting'); ylabel('pearsons R (to F)')
    New_XTickLabel = get(gca,'xtick');
    set(gca,'XTickLabel',New_XTickLabel);

    subplot(2,5,3:4)
    scatter(ssdSort,covFX(ssdSortI),'.')
    xlim([.1 10])
    set(gca,'xscale','log'); grid on
    xlabel('cell variability')
    New_XTickLabel = get(gca,'xtick');
    set(gca,'XTickLabel',New_XTickLabel);
    set(gca,'YTickLabel',[])

    subplot(2,5,5)
    histogram(covFX,'Orientation','horizontal')
    set(gca,'YTickLabel',[])
    grid on

    subplot(2,5,6:7)
    plot(wSort,percInc(wSortI),'.')
    xlim([.1 10])
    set(gca,'xscale','log'); grid on
    xlabel('weighting');  ylabel('activation change (%)')
    New_XTickLabel = get(gca,'xtick');
    set(gca,'XTickLabel',New_XTickLabel);

    subplot(2,5,8:9)
    scatter(ssdSort,percInc(ssdSortI),'.')
    xlim([.1 10])
    set(gca,'xscale','log'); grid on
    xlabel('cell variability')
    New_XTickLabel = get(gca,'xtick');
    set(gca,'XTickLabel',New_XTickLabel);
    set(gca,'YTickLabel',[])

    subplot(2,5,10)
    histogram(percInc,'Orientation','horizontal')
    set(gca,'YTickLabel',[])
    grid on

    set(findall(gcf,'-property','FontSize'),'FontSize',18)   
end

%% function to get 95% confidence interval returns mean m, ci c
function [m,c] = getCI(f)
    % handle cases where we only used one monte carlo sim
    if size(f,1)==1
        m = f;
        c = zeros(size(m));
        return
    end
    % otherwise calculate mean and 95% CIs as usual
    m = mean(f);
    n = size(f,1);
    sem = std(f)/sqrt(n);
    ci95 = tinv(0.025, n-1);
    c = bsxfun(@times, sem, ci95(:));
end