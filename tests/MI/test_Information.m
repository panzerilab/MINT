classdef test_Information < matlab.unittest.TestCase
    methods (Test)

        function testEmptyInput(testCase)
            testCase.verifyError(@() information(), 'information:notEnoughInput');        
        end

        function testInformationValidInput(testCase)
            rng('default');
            dt = 1/500;
            trialendtime = 0.4;
            t_trial = 0:dt:trialendtime;
            nStimuli = 2;
            nTrials = 100;
            % generate response to stimuli
            maxrate = [50 20];
            rate = nan(nStimuli,length(t_trial));
            for i = 1:nStimuli
                signal = normpdf(t_trial,0.2,0.05);
                rate(i,:) =  maxrate(i) * signal / max(signal);
            end
            % generate responses (spike count)
            R = [];
            S = [];
            spike_train = zeros(nStimuli,nTrials, length(t_trial));
            for i = 1:nStimuli
                for j = 1:nTrials
                    spike_train(i,j,:) = poisson_spike_gen(t_trial, rate(i,:), 0);
                    R = [R sum(spike_train(i,j,:))];
                    S = [S i];
                end
            end
            opts.method = 'dr';
            opts.bias = 'naive';
            opts.bin_methodX = 'eqpop';
            opts.verbose = false;
            opts.supressWarnings = true;
            opts.n_binsX = 2;
            outputs = information(R, S, opts, {'I'});
            testCase.verifyEqual( outputs{1}(1), 0.420, 'AbsTol', 0.05);
        end


        function testInformationInvalidOutput(testCase)
            rng('default');
            dt = 1/500;
            trialendtime = 0.4;
            t_trial = 0:dt:trialendtime;
            nStimuli = 2;
            nTrials = 100;
            % generate response to stimuli
            maxrate = [50 20];
            rate = nan(nStimuli,length(t_trial));
            for i = 1:nStimuli
                signal = normpdf(t_trial,0.2,0.05);
                rate(i,:) =  maxrate(i) * signal / max(signal);
            end
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
            opts.method = 'dr';
            opts.bias = 'pt';
            opts.btsp = 100;
            opts.bin_methodX = 'eqpop';
            opts.n_binsX = 2;
            opts.verbose = false;   
            opts.supressWarnings = true;
            testCase.verifyError(@() information(R, S, opts, {'I', 'I'}), 'information:InvalidOutputOption');               
        end

        function testInformationInvalidInput(testCase)
            rng('default');
            dt = 1/500;
            trialendtime = 0.4;
            t_trial = 0:dt:trialendtime;
            nStimuli = 2;
            nTrials = 100;
            % generate response to stimuli
            maxrate = [50 20];
            rate = nan(nStimuli,length(t_trial));
            for i = 1:nStimuli
                signal = normpdf(t_trial,0.2,0.05);
                rate(i,:) =  maxrate(i) * signal / max(signal);
            end
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
            R(1,198) = NaN;
            opts.method = 'dr';
            opts.bias = 'pt';
            opts.btsp = 100;
            opts.bin_methodX = 'eqpop';
            opts.n_binsX = 2;
            opts.verbose = false; 
            opts.supressWarnings = true;            
            testCase.verifyError(@() information(R, S, opts, {'I'}), 'sanitycheck:NaNInput');               
        end
    end
end