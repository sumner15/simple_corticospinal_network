function [Xf, force, fMax, XVideo] = ...
    simulate_CS(X0,w,maxRate,ssd,target,chemotaxis,nDays,dosage,TNP_dose)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [Xf, force, fMax] = simulateSluggish(X0,w,maxRate,ssd)
% 
% inputs:
% X0 - initial neuron activation in 1xN space where N is number or neurons
% w  - weighting for neuron i
% maxRate - max firing rate for neurons (for saturation function)
% ssd - sluggish/speed vector modifier for each neuron
% target - neurons to target. An empty array '[]' is conventional rehab
% chemotaxis - boolean, turn chemotaxis on/off
% *nDays - how many days to test
% *dosage - number of trials on each day (vector)
% TNP_dose - how often to use TNP intervention (1/TNP_dose trials)
% 
% * these values can also be thought of as sum(dosage) as the total number
% of trials to test & nDays as the resolution of the output (we don't want
% to save the result of every trial if we're testing 10s of 1000s). 
% 
% outputs: 
% Xf - final activation vector (neuron firing rates)
% force - force output (1 x nDays vector)
% fMax - maximum theoretical force output
% XVideo - history of activation rates X 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('TNP_dose','var')
   TNP_dose = 5;                %1/n trials gets TNP feedback
end
minRate = 0;                    %minimum firing rate (if not defined)
N = length(X0);                 %number of neurons
force = zeros(1,nDays);         %data out is max force each day
XVideo = cell(1,nDays);         %save out activation for every day
success = false;                %set flag

if isempty(target)       
    fprintf('Simulating (no targeted feedback)...')
    target = 1:N;
else
    fprintf('Simulating targeted feedback (feedback on %i neurons)', length(target))
end

%initial force production 
F0All = getForce(X0,w,1:N);
F0Target = getForce(X0,w,target);
force(1) = F0All;
XVideo{1} = X0;
%max force possible
fMax = getForce(repmat(maxRate,1,N),w,1:N);
fMax = repmat(fMax,1,nDays);
%counting the number of trial types 
nTarget = 0;
nAll = 0; 

for day = 2:nDays    
    % a number of trials for the given day
    for trial = 1:dosage(day)   
                
        % add stochastic noise on every trial (no chemo), or (if chemo) on
        % trials where the previous trial wasn't a success. If the previous
        % trial was a success, we keep the last vi (no update).
        if ~chemotaxis || (chemotaxis && ~success)
            vi = randn(1,N);    %stochastic noise (normal distribution)                       
        end
        success = false;        %reset flag
        Xi = X0+ssd.*vi;        %new firing rate 
        Xi = max(minRate,Xi);   %set saturation limits on x
        Xi = min(maxRate,Xi);   %set saturation limits on x
        
        %simulate force production (1/n of trials are given feedback)
        FiTarget = getForce(Xi,w,target);   %F for targeted neurons only
        FiAll = getForce(Xi,w,1:N);         %F&X for all neurons
        
        % see if the current trial was successful
        %targeted feedback trial
        if length(target)<N && mod(day,TNP_dose)==0 
            if FiTarget>F0Target    
                success = true;
                nTarget = nTarget+1;       
            end
        %no targeted feedback
        else
            if FiAll>F0All
                success = true;
                nAll = nAll+1;
            end
        end                 
        
        % update activation pattern & force if successful
        if success
            F0All = FiAll;            
            F0Target = FiTarget;            
            X0 = Xi;
        end
    end
    
    %save force for the current day to the force vector   
    XVideo{day} = Xi;
    force(day) = getForce(X0,w,1:N);
   
end
Xf = X0; %returns Xf
fprintf('Done.\n');
fprintf('normal success=%i,  target success=%i\n\n', nAll, nTarget)

% calculate the force produced at this activation vector    
function f = getForce(x,w,target)            
    xTarget = zeros(size(x));
    xTarget(target) = x(target);
    wTarget = zeros(size(w));
    wTarget(target) = w(target);    
    f = mean(wTarget.*xTarget);    
end
end

