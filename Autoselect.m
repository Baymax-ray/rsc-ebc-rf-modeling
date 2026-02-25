function Autoselect(resDir)
files = dir(fullfile(resDir, '*.mat'));
if isempty(files)
    error('No .mat files found in %s', resDir);
end
files(strcmp({files.name}, 'results_naive.mat')) = [];
% ---- build group keys by stripping the trailing "_<taskid>.mat" ----
% Expected filename pattern:
%   <receptive_field_type>_<numClosest>closest_<SGE_TASK_ID>.mat
% Group key = <receptive_field_type>_<numClosest>closest_
groupKey = strings(numel(files),1);

for i = 1:numel(files)
    [~, base, ~] = fileparts(files(i).name);
    % remove the last "_<something>" (task id) part
    % e.g. "results_of_dog_1closest_123" -> "results_of_dog_1closest_"
    key = regexprep(base, '_[^_]+$', '_');
    groupKey(i) = string(key);
end

uKeys = unique(groupKey);

fprintf('Found %d groups in %s\n', numel(uKeys), resDir);

for g = 1:numel(uKeys)
    key = uKeys(g);

    idxFiles = find(groupKey == key);
    fprintf('\nProcessing group: %s  (%d files)\n', key, numel(idxFiles));

    % ---- load first file to initialize ----
    firstPath = fullfile(resDir, files(idxFiles(1)).name);
    S = load(firstPath, 'results_to_save');
    if ~isfield(S,'results_to_save')
        warning('Skip %s (no results_to_save)\n', files(idxFiles(1)).name);
        continue
    end
    results_to_save = S.results_to_save;

    nNeurons = numel(results_to_save);
    results_best(nNeurons) = struct('name',[],'sparseness',[],'loglikelihood',[],'AIC',[],'BIC',[],'params',[]);

    for ne = 1:nNeurons
        results_best(ne).name = results_to_save(ne).name;
        results_best(ne).sparseness = results_to_save(ne).sparseness;
        results_best(ne).loglikelihood = results_to_save(ne).loglikelihood;
        results_best(ne).AIC = results_to_save(ne).AIC;
        results_best(ne).BIC = results_to_save(ne).BIC;
        results_best(ne).params = results_to_save(ne).params;
    end

    % ---- update from remaining files ----
    for k = 1:numel(idxFiles)
        fp = fullfile(resDir, files(idxFiles(k)).name);
        S = load(fp, 'results_to_save');
        if ~isfield(S,'results_to_save'), continue; end
        R = S.results_to_save;

        % build a name->index map for this file (faster than find each time)
        namesR = {R.name};

        for ne = 1:nNeurons
            ii = find(strcmp(namesR, results_best(ne).name), 1);
            if isempty(ii), continue; end

            if results_best(ne).loglikelihood < R(ii).loglikelihood
                results_best(ne).sparseness   = R(ii).sparseness;
                results_best(ne).loglikelihood = R(ii).loglikelihood;
                results_best(ne).AIC          = R(ii).AIC;
                results_best(ne).BIC          = R(ii).BIC;
                results_best(ne).params       = R(ii).params;
            end
        end
    end

    outName = sprintf('%sresults_best.mat', key);
    outPath = fullfile(resDir, outName);
    fprintf('Saving: %s\n', outPath);
    save(outPath, 'results_best');
end

disp('All groups done.');
end