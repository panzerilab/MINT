%
%  This source code is part of:
%  NIT - Neuroscience Information Toolbox
%  Copyright (C) 2020  Roberto Maffulli, Miguel Angel Casal Santiago
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
%%% # MI binned methods tutorials
%%% ## Tutorial 1
%%% Tutorial 1 shows how to use NIT to calculate Mutual Information with 1D rate coding

clear all
close all
clc
rng('default')

dt = 1/500;
trialendtime = 0.4;
t_trial = 0:dt:trialendtime;
nStimuli = 2;
nTrials = 100;

cmap = winter(nStimuli);
% generate response to stimuli
maxrate = [50 20];  % peak rate for each stimulus
rate = nan(nStimuli,length(t_trial));
for i = 1:nStimuli
    signal = normpdf(t_trial,0.2,0.05);
    rate(i,:) =  maxrate(i) * signal / max(signal);  %normalize for peak
end

figure()
subplot(2,1,1)
title("Firing rates responses to stimuli")
hold on
for i=1:nStimuli
    plot(t_trial,rate(i,:),'Linewidth',2,...
        'Color', cmap(i,:),...
        'DisplayName',"Response to stimulus " + i);
end
xlabel('t [s]')
ylabel('Firing rate [Hz]')
set(gca,'TickDir','out')
legend()

% generate responses (spike count)
R = [];
S = [];
for i = 1:nStimuli
    for j = 1:nTrials
        spike_train(i,j,:) = poisson_spike_gen(t_trial, rate(i,:), 0);
        R = [R sum(spike_train(i,j,:))];
        S = [S i];
    end
end

subplot(2,1,2)
title("Raster Plot")
hold on
for i = 1:nStimuli
    for j = 1:nTrials
        plot(t_trial(spike_train(i,j,:) == 1), ((i-1)*nTrials+j)*ones(size(t_trial(spike_train(i,j,:) == 1))), ...
            'LineStyle', 'none', 'Marker', '.', 'MarkerSize',3,...
            'MarkerFaceColor',cmap(i,:),'MarkerEdgeColor',cmap(i,:));
    end
end
xlabel('t [s]')
set(gca,'ytick',[])
set(gca,'TickDir','out')

% define options
opts.method = 'dr';
opts.bias = 'qe';
opts.shuff = 100;
opts.bin_methodX = 'eqpop';
opts.n_binsX = 2;
% calculate MI
outputs = information(R, S, opts, {'I'});
I = outputs{1};

%%
%%% ## Tutorial 2
%%% Tutorial 2 shows how to use NIT to calculate Mutual Information with 2D population rate coding and different binning strategies for the two neurons

clear all
close all
clc
rng('default')

dt = 1/500;
trialendtime = 0.4;
t_trial = 0:dt:trialendtime;
nStimuli = 2;
nTrials = 100;
nNeurons = 2;

cmap = figureproperties(nStimuli); %winter(nStimuli);

% generate mean rates
maxrate(1,1) = 10; % neuron 1, stimulus 1
maxrate(1,2) = 50; % neuron 1, stimulus 2
maxrate(2,1) = 20; % neuron 2, stimulus 1
maxrate(2,2) = 30; % neuron 2, stimulus 2

rate = nan(nStimuli,nNeurons,length(t_trial));
for i = 1:nStimuli
    for j = 1:nNeurons
        signal = normpdf(t_trial,0.2,0.05);
        rate(i,j,:) =  maxrate(j,i) * signal / max(signal);  %normalize for peak
    end
end

% plot responses
figure()
for n=1:nNeurons
    for i=1:nStimuli
        subplot(2,nNeurons,n)
        title("Firing rates responses to stimuli for Neuron " + n)
        hold on
        plot(t_trial,squeeze(rate(i,n,:)),'Linewidth',2,...
            'Color', cmap(i,:),...
            'DisplayName',"Response to stimulus " + i);
        xlabel('t [s]')
        ylabel('Firing rate [Hz]')
        set(gca,'TickDir','out')
        legend()
    end
end

% generate responses (spike count)
R = zeros(nNeurons,nStimuli*nTrials);
S = zeros(1,nStimuli*nTrials);
for i = 1:nStimuli
    for j = 1:nTrials
        S((i-1)*nTrials+j) = i;
        for k = 1:nNeurons
            spike_train(i,j,k,:) = poisson_spike_gen(t_trial, squeeze(rate(k,i,:)), 0);
            R(k,(i-1)*nTrials+j) = sum(spike_train(i,j,k,:));
        end
    end
end

% plot raster plots
for n=1:nNeurons
    for i=1:nStimuli
        subplot(2,nNeurons,nNeurons + n)
        title("Spike raster for Neuron " + n)
        hold on
        for j = 1:nTrials
            plot(t_trial(spike_train(i,j,n,:) == 1), ((i-1)*nTrials+j)*ones(size(t_trial(spike_train(i,j,n,:) == 1))), ...
                'LineStyle', 'none', 'Marker', '.', 'MarkerSize',3,...
                'MarkerFaceColor',cmap(i,:),'MarkerEdgeColor',cmap(i,:));
        end
        xlabel('t [s]')
        set(gca,'ytick',[])
        set(gca,'TickDir','out')
    end
end

% define options
opts.method = 'dr';
opts.bias = 'qe';
opts.shuff = 100;
opts.bin_methodX = {'eqspace' 'eqpop'};  % using two different binning methods for the two dimensions: there is no realistic reason to do it for this specific case but it serves as an example of how to use NIT
opts.n_binsX = 2;
% calculate MI
outputs = information(R, S, opts, {'I'});
I = outputs{1};

%%
%%% ## Tutorial 3
%%% Tutorial 3 shows how to use NIT to calculate Mutual Information breakdown with 2D population rate coding in presence of cross-correlated activity between the teo of neurons

clear all
close all
clc
rng('default')

dt = 1/500;
trialendtime = 1;
t_trial = 0:dt:trialendtime;
nStimuli = 2;
nTrials = 500;
nNeurons = 2;

cmap = winter(nStimuli);

% define rate responses of neurons to stimuli
rate(1,1) = 10; % neuron 1, stimulus 1
rate(1,2) = 1; % neuron 1, stimulus 2
rate(2,1) = 10; % neuron 2, stimulus 1
rate(2,2) = 1; % neuron 2, stimulus 2

% cross-correlated response between the two neurons is obtained by
% generating another, stimulus-dependent, spike train with following peak
% rates
ccrate(1) = 10;
ccrate(2) = 15;

% in order to generate a finite cross-correlogram, the spike trains
% corresponding to ccrate are added to the spike trains of the neurons with
% a random lag, sampled from a norm pdf with std lagstd
lagstd = 0.05;              % seconds

% generate responses (spike count)
R = zeros(nNeurons,nStimuli*nTrials);
S = zeros(1,nStimuli*nTrials);

for i = 1:nStimuli
    for j = 1:nTrials
        S((i-1)*nTrials+j) = i;
        
        % generate common spike train to be added both neurons
        ccspike_train = poisson_spike_gen(t_trial, ccrate(i), 0);
        
        % extract random shift
        lag = normrnd(0,lagstd);    % seconds
        lag = round(lag/dt);        % frames
    
        % shift spike train by lag
        shifted_ccspike_train = zeros(size(ccspike_train));
        if lag >= 0
            shifted_ccspike_train(lag+1:end) = ccspike_train(1:end-lag);
        else
            shifted_ccspike_train(1:end+lag) = ccspike_train(1-lag:end);
        end

        for k = 1:nNeurons
            baseline_spike_train = poisson_spike_gen(t_trial,rate(k,i), 0);
            if k == 1
                spike_train(i,j,k,:) =  baseline_spike_train + ccspike_train;
            else
                spike_train(i,j,k,:) =  baseline_spike_train + shifted_ccspike_train;
            end
            R(k,(i-1)*nTrials+j) = sum(spike_train(i,j,k,:));
        end
    end
end

% plot raster plots
for n=1:nNeurons
    for i=1:nStimuli
        subplot(2,nNeurons,n)
        title("Spike raster for Neuron " + n)
        hold on
        for j = 1:nTrials
            plot(t_trial(spike_train(i,j,n,:) ~= 0),...
                ((i-1)*nTrials+j)*ones(size(t_trial(spike_train(i,j,n,:) ~= 0))), ...
                'LineStyle', 'none', 'Marker', '.', 'MarkerSize',3,...
                'MarkerFaceColor',cmap(i,:),'MarkerEdgeColor',cmap(i,:));
        end
        xlabel('t [s]')
        set(gca,'ytick',[])
        set(gca,'TickDir','out')
        
        subplot(2,nNeurons,n+nNeurons)
        title("Spike count for Neuron " + n)
        hold on
        bar(i,mean(sum(squeeze(spike_train(i,:,n,:)),2)),'FaceColor',cmap(i,:))
        ylim([0 20]);
        xticks([1 2]);
        xlabel("Stimuli");
    end
end

% calculate and plot cross-correlogram
for i=1:nStimuli
    for j=1:nTrials
        [cc(i,j,:), lags] = xcorr(squeeze(spike_train(i,j,1,:)),...
            squeeze(spike_train(i,j,2,:)));
    end
end

figure()
hold on
title("Cross-covariogram between the responses of the neurons")
for i=1:nStimuli
    plot(lags*dt,mean(squeeze(cc(i,:,:)),1) - ...
        xcorr(mean(squeeze(spike_train(i,:,1,:))),...
        mean(squeeze(spike_train(i,:,2,:)))),'Linewidth',2,'Color', cmap(i,:),...
        'DisplayName',"Stimulus " + i);
end

% define options
opts.method = 'dr';
opts.bias = 'qe';
opts.shuff = 100;
opts.bin_methodX = 'eqpop';  % using two different binning methods for the two dimensions
opts.n_binsX = 2;
% calculate information breakdown quantities
outputs = information(R, S, opts, {'I', 'ILIN', 'ISS', 'ICI', 'ICD'});
I = outputs{1};
ILIN = outputs{2};
ISS = outputs{3};
ICI = outputs{4};
ICD = outputs{5};

figure()
bar(categorical({'MI', 'ILIN', 'ISS', 'ICI', 'ICD'}), [I(1) ILIN(1) ISS(1) ICI(1) ICD(1)]);
