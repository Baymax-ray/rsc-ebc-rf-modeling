function Closet_loop(receptive_field_type_list,numClosest,SGE_TASK_ID, nworkers)
    rng(SGE_TASK_ID,"twister"); % set the random seed based on the task ID for reproducibility
    %% load the data
    path='./Demo Data'; %path to the folder containing the data
    boundaries_all=Load_boundaries(path); % all patrick's data lives in a 120*120 box, so the boundaries are fixed
    exclude_list={}; %some neurons are excluded because the location is not correct in that session
    files_names= dir(fullfile(path, 'tracking_data', '*.txt')); % the neuron's name is the name of files in both tracking_data and spike_timestamps, using either one is fine
    files_names = {files_names.name}; % convert struct array to cell array of strings
    nNeurons=length(files_names)-length(exclude_list); %number of neurons to be processed
    for rf_idx = 1:length(receptive_field_type_list)
        results_list{rf_idx}=repmat(struct('name',[],'sparseness',[],'loglikelihood',[],'receptive_field_type',[],...
            'AIC',[],'BIC',[],'params',[],'reFireRate',[]), nNeurons, 1);
    end
    rn=0;
    if isempty(gcp('nocreate'))
    parpool('local', nworkers, 'IdleTimeout',Inf);
    end
    for ne = 1: length(files_names)
        name = files_names{ne};
        name=strrep(name,'.txt','');
        %skip if the rootfk.name is in the exclude list
        if ismember(name, exclude_list)
            continue
        end
        rn=rn+1;
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

        %% set up the receptive field and evaluation function
        % receptive_field_type='angle_distance_gaussian';
        angles = 0:3:357; % 3 degree increments

        %% optimization
        nt = numel(timestamp);
        Egorelative_closest = zeros(nt, numClosest, 2);   % use single if possible
        dist_all            = zeros(nt, numel(angles));
        intersectionPts_all = zeros(nt, numel(angles), 2);
        % ─── heavy constants that every worker needs only to READ ──
        anglesConst     = parallel.pool.Constant(angles);
        boundariesConst = parallel.pool.Constant(boundaries_all);

        parfor t = 1:nt %parallel for ray-tracing is a little faster, can save a few seconds for one neuron
            % Small temporaries live inside each worker
            current_position      = [x_position(t), y_position(t)];
            current_head_direction = head_direction(t);
            % retrieve constant data
            angles_t     = anglesConst.Value;
            boundaries_t = Update_boundaries(boundariesConst.Value, t);
            % ray-tracing
            [dist, intersectionPts] = multiRayTrace2D(...
                current_position, current_head_direction + angles_t, boundaries_t);
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
            Egorelative_closest(t,:,:) = tmp;   % one assignment → OK
        end
        clear anglesConst boundariesConst tmp
        % the ga function's syntax is: [x, fval] = ga(fitnessfcn, nvars, A, b, Aeq, beq, lb, ub, nonlcon, options),
        % or [x, fval] = ga(fitnessfcn, nvars, A, b, Aeq, beq, lb, ub, nonlcon, intcon, options)
        % use [] for empty arguments
        % remember to change lb and ub based on neuron's firing pattern
        clear params
        EgoConst = parallel.pool.Constant(Egorelative_closest);
        firedConst = parallel.pool.Constant(fired);

        for rf_idx = 1:length(receptive_field_type_list)
            clear results_to_save
            receptive_field_type = receptive_field_type_list{rf_idx}; % convert cell to string
            disp(['Processing receptive field type: ', receptive_field_type]);
            results_to_save=results_list{rf_idx};
            results_to_save(rn).receptive_field_type=receptive_field_type;
            results_to_save(rn).name=name;
            results_to_save(rn).sparseness=sparseness;

            objFun=@(p) myScoreFunction(p, EgoConst.Value, firedConst.Value,receptive_field_type); % objective function
            switch receptive_field_type
                case 'angle_distance_gaussian'
                    nParams = 4;  % [sigmaAngle, muAngle, sigmaDist, muDist]
                    lb = [5,  0,    1,  0];  % no negative angles, no negative distance
                    ub = [90, 360,  40,   40];  % no angles greater than 360
                    A=[]; b=[];
                case 'angle_gaussian'
                    nParams = 2;  % [sigmaAngle, muAngle]
                    lb = [5,  0];  % no negative angles
                    ub = [90, 360];  % no angles greater than 360
                    A=[]; b=[];
                case 'distance_gaussian'
                    nParams = 2;  % [sigmaDist, muDist]
                    lb = [1, 0];  % no negative distance
                    ub = [40, 40];  % no angles greater than 360
                    A=[]; b=[];
                case 'single_2d_gaussian'
                    nParams = 3;  % [sigma, egocenter_x, egocenter_y]
                    lb = [0.1, -30, -30,];  
                    ub = [30,  30,  30,];
                    A=[]; b=[];
                case 'stretched_2d_gaussian'
                    nParams = 5;  % [sigma_major, ratio, ellipse_angle, egocenter_x, egocenter_y]
                    lb = [1,  0.01,  0,  -30,  -30]; %sigma_minor<sigma_major
                    ub = [ 60,  1,  180,  30,   30];
                    A=[]; b=[]; 
                case 'dog'
                    nParams = 5;  % [ratio, sigma2, a2, egocenter_x, egocenter_y]
                    lb = [0.01, 1, 0, -30, -30]; %sigma1<sigma2
                    ub = [1, 40, 1, 30, 30];  %a2<1
                    A=[]; b=[]; 
                case 'bounded_gaussian'
                    nParams = 4;  % [prefAngle, halfWidth, sigmaDist, muDist]
                    lb = [0,  0,    1,  0];  % no negative angles, no negative distance
                    ub = [360, 90,   40,  40];  % no angles greater than 360
                    A=[]; b=[];
            end
            nParams=nParams+2; %add alpha and offset
            if numClosest<5
                lb = [lb, 0, -6];
                ub = [ub, 5, 6]; 
            elseif numClosest<30
                lb = [lb, 0, -6];
                ub = [ub, 3, 6];
            else
                lb = [lb, 0, -6];
                ub = [ub, 1, 12];
            end
            if ~isempty(A), A = [A, zeros(size(A, 1), 2)]; end
                %% ── hybrid local-search options (fmincon) ──────────────────────────────
            hybridopts = optimoptions('fmincon', ...
            'Algorithm',               'interior-point', ...  % robust default
            'Display',                 'none',        ...     % quiet; change to 'iter' if debugging
            'FiniteDifferenceType',    'central',     ...     % better gradient estimate
            'FiniteDifferenceStepSize',1e-4,         ...     % small step – avoids angle seam
            'MaxIterations',           5e3,           ...
            'MaxFunctionEvaluations',  5e4);
            popSize   = max( 400, 40*nParams );   % ≈40×variables is a good heuristic
            eliteCnt  = ceil( 0.05 * popSize );   % preserve top 5 %
            options = optimoptions('ga', ...
                'UseParallel',          true, ...
                'PopulationSize',       popSize, ...
                'EliteCount',           eliteCnt, ...
                'CrossoverFraction',    0.8, ...
                'MutationFcn',          {@mutationadaptfeasible}, ...
                'InitialPopulationRange',[lb; ub], ...   % helps GA start inside bounds
                'MaxGenerations',       150,     ...
                'FunctionTolerance',    1e-6,    ...
                'MaxStallGenerations',  40,      ...
                'Display',              'off',   ...     % GA console output
                'HybridFcn',            {@fmincon, hybridopts}, ...
                'OutputFcn',            @myOutputFcn);
            % Run the genetic algorithm
            [p, fval] = ga(objFun, nParams, A, b, [], [], lb, ub, [], options);
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
            %display the optimized parameters
            par_print=sprintf('%.3f ', p);
            fprintf('Optimized parameters: %s\n', par_print);
            disp("Best score: ");
            loglikelihood= -fval;
            disp(loglikelihood);
            results_to_save(rn).loglikelihood=loglikelihood;
            results_to_save(rn).params=p;
            fprintf('log-likelihood of the simulated firing: %.3f\n', loglikelihood);
            AIC=2*nParams-2*loglikelihood; %Akaike Information Criterion (lower is better)
            BIC=log(length(timestamp))*nParams-2*loglikelihood; %Bayesian Information Criterion (lower is better)
            fprintf('AIC: %.3f\n', AIC);
            fprintf('BIC: %.3f\n', BIC);
            results_to_save(rn).AIC=AIC;
            results_to_save(rn).BIC=BIC;
            results_list{rf_idx} = results_to_save;
        end
    end
    %% save the results
    SGE_TASK_ID = num2str(SGE_TASK_ID); % ensure SGE_TASK_ID is a string
    for rf_idx = 1:length(receptive_field_type_list)
        receptive_field_type = receptive_field_type_list{rf_idx}; % convert cell to string
        results_to_save=results_list{rf_idx};
        disp([ 'Saving the results of ' receptive_field_type ' to results_of ' receptive_field_type SGE_TASK_ID '_closest.mat']);
        save(['./Demo results/results_of_' receptive_field_type '_' num2str(numClosest) 'closest_' SGE_TASK_ID '.mat' ],'results_to_save')
        disp('Done!');        
    end
end

%% functions
function score=myScoreFunction(p, Egorelative_closest, fired, receptive_field_type)
    switch receptive_field_type
        case 'angle_distance_gaussian'
            param.sigmaAngle = p(1);
            param.muAngle    = p(2);
            param.sigmaDist  = p(3);
            param.muDist     = p(4);
        case 'angle_gaussian'
            param.sigmaAngle = p(1);
            param.muAngle    = p(2);
        case 'distance_gaussian'
            param.sigmaDist  = p(1);
            param.muDist     = p(2);
        case 'single_2d_gaussian'
            param.sigma = p(1);
            param.egocenter = [p(2), p(3)];
        case 'stretched_2d_gaussian'
            param.sigma_major = p(1);
            param.ratio = p(2);
            param.ellipse_angle = p(3);
            param.egocenter = [p(4), p(5)];
        case 'dog'
            param.ratio = p(1);
            param.sigma2 = p(2);
            param.a2 = p(3);
            param.egocenter = [p(4), p(5)];
        case 'bounded_gaussian'
            param.prefAngle = p(1);
            param.halfWidth = p(2);
            param.sigmaDist = p(3);
            param.muDist    = p(4);
    end
    param.alpha=p(end-1);
    param.offset=p(end);
    [N, X, ~] = size(Egorelative_closest); %X is the number of closest points
    pts = reshape(Egorelative_closest, [N*X, 2]);
    rfPerPoint=computeTotalReceptiveField(receptive_field_type, pts, param); %calculate the raw receptive field value
    rfPerPoint = reshape(rfPerPoint, [N, X]); %reshape it back to N*X
    rfValue = sum(rfPerPoint, 2); %sum the receptive field value of the closest points
    %activation function
    weight= param.alpha.*rfValue-param.offset; %calculate the weight
    % Calculate the score as the negative log-likelihood
    p = 1 ./ (1 + exp(-weight)); % logistic function
    % avoid exact 0 or 1
    epsval = 1e-15;  
    p = min( max(p, epsval), 1 - epsval );
    score=-sum(fired.*log(p)+(1-fired).*log(1-p)); %negative log-likelihood
end


function [state, options, optchanged] = myOutputFcn(options, state, ~)
    % Find best individual in the current generation
    [bestScore, bestIdx] = min(state.Score);
    bestParam = state.Population(bestIdx, :);

    % Only print every 10 generations
    if mod(state.Generation, 10) == 0
        % Print generation #
        fprintf('Generation %d: ', state.Generation)

        % Print parameters with 4 decimal places
        fprintf('Best params = [');
        fprintf('%.4f ', bestParam);    % change %.4f to %.2f, %.6f, etc. as you like
        fprintf(']  ');

        % Print the score
        fprintf('Score = %.4f\n', bestScore);
    end

    optchanged = false;
end