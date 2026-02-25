function simulate_draw(monte_carlo_simulations,draw,fps)
% It is used to simulate the drawing of receptive fields based on the _best.mat files
folder = './Demo results'; % folder containing the _best.mat files
files = dir(fullfile(folder, '*results_best.mat'));  % get all results_best.mat files


if monte_carlo_simulations==0 && draw==0
    warning('Both monte_carlo_simulations and draw are set to 0, nothing will be done.');
    return;
end

if isempty(files)
        fprintf('No *_best.mat files found in: %s\n', folder);
        return;
end

for k = 1:numel(files)
    fname = files(k).name;
    fpath = fullfile(folder, fname);
    tokens = regexp(fname, '^results_of_(?<type>.+?)_(?<num>\d+)closest_.*_best\.mat$', 'names');
    time = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');
    fprintf('Processing file: %s at %s\n', fname, char(time));
    if isempty(tokens)
        warning('Skipping (name not recognized): %s', fname);
        continue;
    end
    receptive_field_type = tokens.type;
    pts_number = str2double(tokens.num);
    S = load(fpath);
    R=S.results_best;
    results_best=Simulate_helper(R, receptive_field_type, pts_number, draw, monte_carlo_simulations, fps);
    newName = regexprep(fname, 'results_best\.mat$', 'results_best_filled.mat');
    newPath = fullfile(folder, newName);
    time = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');
    fprintf('Finished processing at %s\n', char(time));
    if monte_carlo_simulations~=0
        save(newPath, 'results_best');
        fprintf('Filled: %s\n', newName);
    end
end
end

function S_out=Simulate_helper(Result, receptive_field_type, numClosest, draw, monte_carlo_simulations, fps)
    rng(12345,"twister"); 
    %% load the data
    path='./Demo Data'; %path to the folder containing the data
    boundaries_all=Load_boundaries(path); % all patrick's data lives in a 120*120 box, so the boundaries are fixed
    exclude_list={}; %some neurons are excluded because the location is not correct in that session
    files_names= dir(fullfile(path, 'tracking_data', '*.txt')); % the neuron's name is the name of files in both tracking_data and spike_timestamps, using either one is fine
    files_names = {files_names.name}; % convert struct array to cell array of strings
    nNeurons=length(files_names)-length(exclude_list); %number of neurons to be processed
    if nNeurons~=length(Result)
        error('Number of neurons does not match: %d (data) vs %d (model)', nNeurons, length(Result));
    end
    S_out(nNeurons)=struct('name',[],'sparseness',[],'loglikelihood',[],'pvalue',[],'AIC',[],'BIC',[],'params',[],'reFireRate',[]);
    rn=0;
    for ne = 1 : length(files_names)
        name = files_names{ne};
        name=strrep(name,'.txt','');
        nidx = find(strcmp({Result.name}, name), 1);
        if isempty(nidx)
            continue
        end
        rn=rn+1;
        S_out(rn).name=name;
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
        if sparseness ~= Result(nidx).sparseness
            warning('Sparseness does not match for %s: %.4f (data) vs %.4f (model)', name, sparseness, Result(nidx).sparseness);
        end
        S_out(rn).sparseness=sparseness;

        angles = 0:3:357; % 3 degree increments
        %% ray tracing
        nt = numel(timestamp);
        Egorelative_closest = zeros(nt, numClosest, 2);   % use single if possible
        dist_all            = zeros(nt, numel(angles));
        intersectionPts_all = zeros(nt, numel(angles), 2);
        for t = 1:nt %parallel for ray-tracing is a little faster, can save a few seconds for one neuron
            % Small temporaries live inside each worker
            current_position      = [x_position(t), y_position(t)];
            current_head_direction = head_direction(t);
            % retrieve constant data
            boundaries_t = Update_boundaries(boundaries_all, t);
            % ray-tracing
            [dist, intersectionPts] = multiRayTrace2D(...
                current_position, current_head_direction + angles, boundaries_t);
            % write whole row-slices ⇒ "sliced variables" are legal
            dist_all(t , :)      = dist;
            intersectionPts_all(t, :, :) = intersectionPts;
            % egocentric conversion for the K closest points
            [~, sorted_index] = sort(dist);
            R = [ cosd(current_head_direction)  sind(current_head_direction) ;
                -sind(current_head_direction)  cosd(current_head_direction) ];
            tmp = zeros(numClosest,2);
            for k = 1:numClosest
                idx      = sorted_index(k);
                closestP = intersectionPts(idx, :);
                tmp(k,:) = (R * (closestP' - current_position'))';
            end
            Egorelative_closest(t,:,:) = tmp;
        end

        %% recover params
        p=Result(nidx).params;
        S_out(rn).loglikelihood=Result(nidx).loglikelihood;
        S_out(rn).AIC=Result(nidx).AIC;
        S_out(rn).BIC=Result(nidx).BIC;
        S_out(rn).params=p;
        switch receptive_field_type
            case 'angle_distance_gaussian'
                params.sigmaAngle = p(1);
                params.muAngle    = p(2);
                params.sigmaDist  = p(3);
                params.muDist     = p(4);
            case 'angle_gaussian'
                params.sigmaAngle = p(1);
                params.muAngle    = p(2);
            case 'distance_gaussian'
                params.sigmaDist  = p(1);
                params.muDist     = p(2);
            case 'single_2d_gaussian'
                params.sigma = p(1);
                params.egocenter = [p(2), p(3)];
            case 'stretched_2d_gaussian'
                params.sigma_major = p(1);
                params.ratio = p(2);
                params.ellipse_angle = p(3);
                params.egocenter = [p(4), p(5)];
            case 'dog'
                params.ratio = p(1);
                params.sigma2 = p(2);
                params.a2 = p(3);
                params.egocenter = [p(4), p(5)];
            case 'bounded_gaussian'
                params.prefAngle = p(1);
                params.halfWidth = p(2);
                params.sigmaDist = p(3);
                params.muDist    = p(4);
        end
        params.alpha=p(end-1);
        params.offset=p(end);
        %% visualize the receptive field
        if draw ~=0
            xMin = -30; xMax = 30;  nx = 500;
            yMin = -30; yMax = 30;  ny = 500;
            xRange = linspace(xMin, xMax, nx);
            yRange = linspace(yMin, yMax, ny);
            [X, Y] = meshgrid(xRange, yRange);
            RF = zeros(size(X));
            fake_head_direction=90; %head to the top
            fake_position=[0,0]; %center at the origin
            for i = 1:nx
                for j = 1:ny
                    target_position = [X(i,j), Y(i,j)];
                    w = params.alpha*Check_receptive_field(receptive_field_type, fake_position, fake_head_direction, target_position, params)-params.offset;
                    RF(i,j) = 1/(1+exp(-w)); % logistic function
                end
            end
            f = figure('Visible','off');
            % Use imagesc to display the RF values. 'axis xy' ensures that y increases upward.
            imagesc(xRange, yRange, RF);
            axis xy;  
            colormap(jet);  
            colorbar;
            xlabel('Right (+) / Left (-)  [cm]');
            ylabel('Forward (+) / Backward (-)  [cm]');
            title('Receptive Field (Egocentric)');
            hold on;
            plot(fake_position(1), fake_position(2), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'black'); % plot the position of the animal
            quiver(fake_position(1), fake_position(2), 10*cosd(fake_head_direction), 10*sind(fake_head_direction), 'black', 'LineWidth', 2,'MaxHeadSize', 2); % plot the head direction of the animal
            hold off;
            % convert numClosest to string for the title
            folderPath = fullfile('.', 'Demo figures', [num2str(numClosest) 'pts'], receptive_field_type);
            if ~exist(folderPath, 'dir')
                mkdir(folderPath);
            end

            fileName = [name '_' receptive_field_type '_' num2str(numClosest) 'closest_rf.png'];
            fullPath = fullfile(folderPath, fileName);
            saveas(f, fullPath);
            close(f); % close the figure to save memory
        end
        %% simulate with parameters
        if monte_carlo_simulations~=0
            prob_simulated=zeros(length(timestamp),1);
            closest_dist=zeros(length(timestamp),1);
            closest_angle=zeros(length(timestamp),1);
            loglikelihood=zeros(length(timestamp),1);
            for t=1:length(timestamp)
                %get the current position and head direction
                current_position=[x_position(t),y_position(t)];
                current_head_direction=head_direction(t);
                %check the receptive field to see if it should fire
                dist = dist_all(t,:); %distance to the closest point
                intersectionPts = squeeze(intersectionPts_all(t,:,:)); %closest point on the boundary
                [~, sorted_index] = sort(dist);
                topIdx = sorted_index(1:numClosest);
                closest_Xpts = intersectionPts(topIdx, :);
                closest_dist(t)=dist(sorted_index(1));
                closest_angle(t)=angles(sorted_index(1));
                weight=params.alpha * Check_receptive_field(receptive_field_type, current_position, current_head_direction, closest_Xpts, params) - params.offset; %calculate the weight
                prob=1/(1+exp(-weight)); % logistic function
                % avoid exact 0 or 1
                epsval = 1e-15;  
                prob = min( max(prob, epsval), 1 - epsval );
                loglikelihood(t)=fired(t)*log(prob)+(1-fired(t))*log(1-prob);
                prob_simulated(t)=prob;
            end
        %% plot the simulated firing
            small_size=1;
            large_size=20;
            prob=ones(size(timestamp))*sparseness; %probability of firing=average firing rate
            %randomly generate the firing based on the probability of firing
            fired_simulated=rand(size(prob_simulated))<prob_simulated;
            %permutation test
            loglikelihood_random=zeros(2000,1);
            for i=1:2000
                %randomly permute the firing
                fired_random=rand(size(timestamp))<prob;
                %calculate the loglikelihood of the random firing
                loglikelihood_random(i)=sum(fired_random.*log(prob_simulated)+(1-fired_random).*log(1-prob_simulated));
            end
            pvalue=length(find(loglikelihood_random>=sum(loglikelihood)))/2000; %pvalue of the loglikelihood
            S_out(rn).pvalue=pvalue;
            re_ave_fire=mean(fired_simulated)*fps; %regenerated average firing rate in Hz
            S_out(rn).reFireRate=re_ave_fire;
        end
        if draw ~=0 && monte_carlo_simulations~=0
        % plot the Regenerated firing
            f=figure('Visible','off');
            scatter(x_position, y_position, small_size+fired_simulated*large_size, head_direction, 'filled');
            colormap(hsv);           % Suitable for directional data (0 to 360)
            clim([0 360]);
            xlim([0 120]);
            ylim([0 120]);
            xlabel('X Position');
            ylabel('Y Position');
            title( sprintf('Regenerated Firing (fire rate=%.3f Hz, p=%.3f)', re_ave_fire, pvalue) );
            % Add the directional legend
            Add_directional_legend([0.9 0.885 0.09 0.09]);
            saveas(f, ['./Demo figures/' num2str(numClosest) 'pts/' receptive_field_type '/' name '_' receptive_field_type '_' num2str(numClosest) 'closest_regenerated_firing.png']);
            close(f);
        end
    end
end