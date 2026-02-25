function Basics()
rng(12345,"twister");
%% load the data
path='./Demo Data'; %path to the folder containing the data
boundaries_all=Load_boundaries(path); % all patrick's data lives in a 120*120 box, so the boundaries are fixed
exclude_list={}; %some neurons are excluded because the location is not correct in that session
files_names= dir(fullfile(path, 'tracking_data', '*.txt')); % the neuron's name is the name of files in both tracking_data and spike_timestamps, using either one is fine
files_names = {files_names.name}; % convert struct array to cell array of strings
nNeurons=length(files_names)-length(exclude_list); %number of neurons to be processed
results_to_save(nNeurons)= struct('name',[],'sparseness',[],'loglikelihood',[],'pvalue',[],'AIC',[],'BIC',[],'params',[],'reFireRate',[]);
rn=0;
fps=30; %frame per second of the tracking data
%% process each neuron
for ne = 1: length(files_names)
    name = files_names{ne};
    name=strrep(name,'.txt','');
    %skip if the rootfk.name is in the exclude list
    if ismember(name, exclude_list)
        continue
    end
    rn=rn+1;
    results_to_save(rn).name=name;
    fprintf('\nProcessing %s\n', name);
    [timestamp,x_position,y_position,head_direction]=Load_tracking_data(path,['tracking_data/' name '.txt']);
    spike_times_real=Load_spike_times(path,['spike_timestamps/' name '.txt']);
    fired = false(size(timestamp)); %if this frame is fired
    fired_count=zeros(size(timestamp)); %how many times this frame is fired
    for i = 1:length(spike_times_real) % find the closest timestamp to the real spike time
        [~, idx] = min(abs(timestamp - spike_times_real(i)));
        fired(idx) = true;
        fired_count(idx)=fired_count(idx)+1;
    end
    %use the min and max of x and y to define the boundaries of the environment because that data is not recorded in the initial file
    %if any of x, y, or head direction is NaN or < 0 or > 120, remove that index from all of them
    nanIndex = isnan(x_position) | isnan(y_position) | isnan(head_direction) | x_position < 0 | x_position > 120 | y_position < 0 | y_position > 120;
    x_position(nanIndex) = [];
    y_position(nanIndex) = [];
    head_direction(nanIndex) = [];
    fired(nanIndex) = [];
    fired_count(nanIndex) = [];
    timestamp(nanIndex) = [];
    head_direction = mod(head_direction,360); % head direction should be in [0, 360]
    sparseness=sum(fired)/length(fired);
    fprintf('sparseness: %.6f\n', sparseness);
    results_to_save(rn).sparseness=sparseness;
    
    %% plot the raw firing
    small_size=1;
    large_size=20;
    f = figure('Visible','off');
    marker_size=small_size*ones(size(timestamp));
    marker_size(fired)=small_size+large_size;
    scatter(x_position, y_position, marker_size, mod(head_direction,360), 'filled');
    colormap(hsv);           % Suitable for directional data (0 to 360)
    clim([0 360]); 
    xlim([0 120]);
    ylim([0 120]);
    xlabel('X Position');
    ylabel('Y Position');
    title(sprintf('Raw Firings (fire rate=%.3f Hz)', sparseness*fps));
    % Add the directional legend
    Add_directional_legend([0.9 0.885 0.09 0.09]);
    saveas(f, ['.\Demo figures\' name '_raw_firing.png']);
    close(f);
    

    %% fit the most simple model as a baseline
    prob=ones(size(timestamp))*sparseness; %probability of firing=average firing rate
    loglikelihood=sum(fired.*log(prob)+(1-fired).*log(1-prob)); %loglikelihood of the model
    AIC=-2*loglikelihood+2*1; %AIC=-2*loglikelihood+2*k, k is the number of parameters
    BIC=-2*loglikelihood+log(length(timestamp))*1; %BIC=-2*loglikelihood+log(n)*k, n is the number of data points
    fprintf('AIC: %.6f\n', AIC);
    fprintf('BIC: %.6f\n', BIC);
    results_to_save(rn).loglikelihood=loglikelihood;
    results_to_save(rn).AIC=AIC;
    results_to_save(rn).BIC=BIC;
    %regenerate the fire plots
    f = figure('Visible','off');
    fired_simulated=rand(size(timestamp))<prob; %simulate the firing based on the probability
    re_ave_fire=mean(fired_simulated)*fps; %regenerated average firing rate in Hz
    results_to_save(rn).reFireRate=re_ave_fire;
    fprintf('Re-ave fire: %.6f\n', re_ave_fire);
    marker_size=small_size*ones(size(timestamp));
    marker_size(fired_simulated)=small_size+large_size;
    scatter(x_position, y_position, marker_size, mod(head_direction,360), 'filled');
    colormap(hsv);           % Suitable for directional data (0 to 360)
    clim([0 360]); 
    xlim([0 120]);
    ylim([0 120]);
    xlabel('X Position');
    ylabel('Y Position');
    title( sprintf('Regenerated Firing (fire rate=%.3f Hz)', re_ave_fire) );
    % Add the directional legend
    Add_directional_legend([0.9 0.885 0.09 0.09]);
    saveas(f, ['.\Demo figures\' name '_naive_regenerated_firing.png']);
    close(f);
end
%% save the results to a mat file
disp([ 'Saving the results of all neurons to ' 'results_naive.mat'])
save('.\Demo results\results_naive.mat','results_to_save')
disp('Done!')
end

