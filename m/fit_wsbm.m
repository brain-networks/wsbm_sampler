clear all
close all
clc

addpath(genpath('WSBM_v1.2'));      % add wsbm toolbox to path (http://tuvalu.santafe.edu/~aaronc/wsbm/)
load ../data/A                      % load connectivity matrix
N = length(A);                      % number of nodes
A(1:(N + 1):end) = 0;               % remove diagonal

W_distr = 'normal';                 % distribution over edge weights
alph = 0;                           % balance weight/edge probability distributions (alph = 0 means we conly care about weights)
muMaxIter = 250;                    % number of optimization steps for mu parameter
mainMaxIter = 250;                  % number of optimization steps for main loop
mainTol = 1.0000e-03;               % tolerance -- smaller = better fit, but longer runtime
muTol = 1.0000e-03;                 % tolerance -- smaller = better fit, but longer runtime

k = 6;                              % number of communities to detect
edgeList = Adj2Edg(A);              % convert matrix to edges
mx = inf;                           % initialize number of communities to inf
while mx ~= k                       % run algorithm until you get k communities -- due to "tie-breaking" you can get fewer than k
    
    tic
    [ci,m] = wsbm(edgeList,k,...    % run model
        'NumTrials',1,...
        'muMaxIter',muMaxIter,...
        'mainMaxIter',mainMaxIter,...
        'W_distr',W_distr,...
        'mainTol',mainTol,...
        'muTol',muTol,...
        'alpha',alph);
    mx = max(length(unique(ci)));   % number of unique communities
    toc
end
[C,Morph] = fcn_comm_motifs(A,ci);  % get community motif types -- type 'help fcn_comm_motifs' in command line for more information