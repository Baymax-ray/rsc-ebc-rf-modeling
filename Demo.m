%%
% This is a demo file for the article "Computational Structure of Egocentric Boundary Cell Responses in Retrosplenial Cortex"
% Run this demo will reproduce Figure 1 in the paper.

% rf is the list of receptive field models that we used in the paper.
%rf = {'angle_distance_gaussian','angle_gaussian','distance_gaussian','single_2d_gaussian','stretched_2d_gaussian','dog','bounded_gaussian'}; 
rf = {'angle_distance_gaussian','dog','bounded_gaussian'}; 
% bps is the amount of boundary points included, 1 means we only use the closest boundary point for this demo, could be changed to higher values to include more boundary points in the model fitting. (max 120 for our script, which includes all boundary points visible in the arena)
bps = 1 ; 
% r is the amount of repetitions for the model fitting. It was set to 5 in the paper. The higher the value, the more stable the model fitting results, but also the longer the running time.
r =5;
% nworkers is the number of parallel workers for the parpool. The higher the value, the faster the running time, but also the more computational resources needed. You can adjust it based on your own computational resources.
nworkers = 4;

monte_carlo_simulations=1; % whether to do monte carlo simulations for p-value calculation
draw=1; % 1 to really draw the simulation, 0 only simulate without storing the figures
fps=30; %frame per second of the tracking data, used for calculating the firing rate in the figures. It does not affect the model fitting results.


% fiting results will be saved in the folder './Demo results'. You can change the path if you want to save the results in a different location.
outDir = './Demo results';
if ~exist(outDir,'dir')
    mkdir(outDir);
end
% figures will be saved in the folder './Demo figures'. You can change the path if you want to save the figures in a different location.
figDir = './Demo figures';
if ~exist(figDir,'dir')
    mkdir(figDir);
end
%figures will be sorted into different subfolders based on the amount of boundary points included in the model fitting.
figDir2 = ['./Demo figures/',num2str(bps),'pts'];
if ~exist(figDir2,'dir')
    mkdir(figDir2);
end

%% run the model fitting in parallel for different repetitions and different receptive field models. The results will be saved in the folder './Demo results' with the name format 'receptive_field_type_numClosest_closest_repetitionID.mat'.
for i=1:r
    Closet_loop(rf,bps,i,nworkers);
end

% Negtive receptive field models are fitted in a similar way, but allowing negative parameters.
% for i=r+1:2*r
%     Closet_loop_Neg(rf,bps,i,nworkers);
% end



%% Auto-select best fitting for ALL result groups in ./results/ and save as results_best.mat in the same folder.
Autoselect(outDir);

%% figure generation

Basics(); % raw firing and random firing
Simulate_draw(monte_carlo_simulations,draw,fps); % simulated firing based on the fitted model parameters