function useParallel = configureDemoParallelPool(nWorkers)
%CONFIGUREDEMOPARALLELPOOL Open or reuse a Windows local parallel pool.

useParallel = false;
nWorkers = max(1, round(nWorkers));
if nWorkers <= 1
    fprintf('Parallel execution disabled (nWorkers = 1).\n');
    return
end

hasParallel = license('test', 'Distrib_Computing_Toolbox') && ...
    exist('parpool', 'file') == 2 && exist('gcp', 'file') == 2;
if ~hasParallel
    warning('Parallel Computing Toolbox is unavailable; using serial loops.');
    return
end

try
    pool = gcp('nocreate');
    if isempty(pool)
        pool = parpool('local', nWorkers);
    elseif pool.NumWorkers ~= nWorkers
        warning('Using the existing pool with %d workers (requested %d).', ...
            pool.NumWorkers, nWorkers);
    end
    useParallel = ~isempty(pool);
catch ME
    warning('Could not open a local parallel pool; using serial loops: %s', ...
        ME.message);
end
end
