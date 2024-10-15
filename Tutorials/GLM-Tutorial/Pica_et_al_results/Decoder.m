clc, clear,
A1_B67_0630 = load('A1_B67_0630_exp_data.mat');

A1_B67_0630 = A1_B67_0630.exp_data;

predictors = struct();

%% Uncoupled model predictors
% Only task-predictors are used
% - Velocity, Position, Stimuli, ViewAngle, Choice, Time for reward/error

for trial = 1:length(A1_B67_0630)
    stimulus_onFrame = A1_B67_0630(trial).frame_timing_info.stimulus_on_frames;
    stimulus_offFrame = A1_B67_0630(trial).frame_timing_info.stimulus_end_frames;
    reward_onFrame = A1_B67_0630(1).frame_timing_info.reward_on_frames;
    reward_offFrame = A1_B67_0630(1).frame_timing_info.reward_end_frames;

    nFrames = length(A1_B67_0630(trial).y_position);

    A1_B67_0630(trial).stimulus_onset = zeros(1, nFrames);
    A1_B67_0630(trial).reward_onset = zeros(1, nFrames);

    for interval = 1:numel(stimulus_onFrame)
        A1_B67_0630(trial).stimulus_onset(stimulus_onFrame(interval):stimulus_offFrame(interval)) = 1;
    end

    A1_B67_0630(trial).reward_onset(reward_onFrame:reward_offFrame) = 1;

    time = (0:1:nFrames-1);

    % Stimulus onset
    numBasisFunctions_sound = 12;
    FWHM_sound = (170/1000)*15.6;
    sigma_sound = FWHM_sound / (2 * sqrt(2 * log(2)));
    duration_sound = 2 * 15.6;
    numRepeats_sound = 3;
    numSoundLocations = 8;
    allBasisFunctions_sound = zeros(numBasisFunctions_sound * numRepeats_sound * numSoundLocations, nFrames);

    for direction = 1:numSoundLocations

        for repeat = 1:numRepeats_sound

            index_sound = ((direction - 1) * numRepeats_sound + (repeat - 1)) * numBasisFunctions_sound + 1;
            
            if numel(stimulus_onFrame) >= repeat && direction == A1_B67_0630(trial).sound_location

                sound_onset = stimulus_onFrame(repeat);
                timePoints_sound = linspace(sound_onset, sound_onset+duration_sound, numBasisFunctions_sound);


                for i = 1:numBasisFunctions_sound
                     mu_sound = round(timePoints_sound(i));
                     gauss_function = exp(-(((time - mu_sound).^2) / (2 * sigma_sound^2)));
                     allBasisFunctions_sound(index_sound + i - 1, :) = gauss_function;
                end

            else
                for i = 1:numBasisFunctions_sound
                    allBasisFunctions_sound(index_sound + i - 1, :) = zeros(1, numel(time));
                end
            end
        end
    end
    predictors(trial).BasisFunctions_sound = allBasisFunctions_sound;

    % Reward/Error
    FWHM_ce = (500/1000)*15.6;
    sigma_ce = FWHM_ce / (2 * sqrt(2 * log(2)));
    duration_ce = 2 * 15.6;
    numBasisFunctions_ce = 8;
    allBasisFunctions_ce = zeros(numBasisFunctions_ce,  nFrames);
    timePoints_ce = linspace(reward_onFrame, reward_onFrame+duration_ce, numBasisFunctions_ce);

    if A1_B67_0630(trial).correct == 1
        for i = 1:4
            mu_ce = round(timePoints_ce(i));
            gauss_function = exp(-(((time - mu_ce).^2) / (2 * sigma_ce^2)));
            allBasisFunctions_ce(i, :) = gauss_function;
        end
        for i = 5:8
            allBasisFunctions_ce(i, :) = zeros(1, numel(time));
        end 
    else
        for i = 1:4
             allBasisFunctions_ce(i, :) = zeros(1, numel(time));
        end
        for i = 5:8
            mu_ce = round(timePoints_ce(i));
            gauss_function = exp(-(((time - mu_ce).^2) / (2 * sigma_ce^2)));
            allBasisFunctions_ce(i, :) = gauss_function;
        end 
    end
    predictors(trial).BasisFunctions_ce = allBasisFunctions_ce;

    % Running velocity
    FWHM_rv = (240/1000)*15.6;
    sigma_rv = FWHM_rv / (2 * sqrt(2 * log(2)));
    duration_rv = 1 * 15.6;
    numBasisFunctions_rv = 8;
    numRunningDirections = 4;

    allBasisFunctions_rv = zeros(numBasisFunctions_rv * numRunningDirections, nFrames);

    for direction = 1:numRunningDirections
        index_speed = (direction - 1) * numBasisFunctions + 1;
        timePoints_speed = linspace(-duration_speed, duration_speed, numBasisFunctions);

        for i = 1:numBasisFunctions
            mu = timePoints_speed(i);
            gauss_function_speed = exp(-((time - mu).^2) / (2 * sigma^2));
            allBasisFunctions_speed(index_speed + i - 1, :) = gauss_function_speed;
        end
    end


end


figure
plot(allBasisFunctions_sound);

sigma_reward_error = round((500 / sqrt(2 * log(2)) / 1000)*15.6);
duration_error = round(2 * 15.6);
allBasisFunctions_reward_error = zeros(8, nFrames);

for i = 1:numel(reward_times)
    index_reward_error = (i - 1) * 8 + 1;
    timePoints_reward_error = linspace(0, duration_reward_error, 8);
    allBasisFunctions_reward_error(index_reward_error:(index_reward_error + 7), :) = exp(-(time - reward_times(i) - timePoints_reward_error).^2 / (2 * sigma_reward_error^2));
    allBasisFunctions_reward_error(index_reward_error:(index_reward_error + 7), :) = exp(-(time - error_times(i) - timePoints_reward_error).^2 / (2 * sigma_reward_error^2));
end

function f = gaussian_distribution(t, mu, sigma)
    p = -(1/2) * ((x - mu)/sigma) .^ 2;
    A = 1/(sigma * sqrt(2*pi));
    f = A.*exp(p);
end 

