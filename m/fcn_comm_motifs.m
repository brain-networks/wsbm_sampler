function [C,Morph] = fcn_comm_motifs(A,Ci,verbose)
% fcn_comm_motifs       calculate community motifs
%
%   [C,Morph] = fcn_comm_motifs(A,Ci) takes as input [N x N] connectivity
%   matrix, A, and an [N x NReps] matrix of community labels. This code
%   identifies the motif formed by pairs of communities and maps those
%   motifs back to individual nodes. There are four motif classes:
%   Assorative, Disassortative, Core, and Periphery. Based on a node's
%   participation in different classes, we also calcualte a Diversity index
%   (an entropy) whose value is greatest when a node participates in all
%   motif types uniformly.
%
%   In addition, this code outputs coordinates and labels for a community
%   motif morphospace. Each point represents a pair of communities and is
%   labeled according to their motif type.
%
%   Inputs:
%           
%       A,      connectivity matrix (can be weighted, directed)
%       Ci,     community labels
%
%   Outputs:
%
%       C,      [N x 5] matrix whose columns are node-level measures of
%               assortativity, disassortativity, coreness, periphery-ness,
%               and the diversity index.
%       Morph,  Coordinates and labels for a morphospace
%
%   Example:
%
%   >> % set some parameters and initialize matrix for labels
%   >> NReps = 10; K = 5; N = length(A); edgeList = Adj2Edg(A);
%   >> Ci = zeros(N,NReps);
%   >> % run the WSBM -- you should specify parameters
%   >> for iRep = 1:NReps; Ci(:,iRep) = wsbm(edgeList,K); end;
%   >> % get community motifs
%   >> [C,Morph] = fcn_comm_motifs(A,Ci);
%   >> % draw morphospace -- you could take log of first three coordinates
%   >> scatter3(Morph(:,1),Morph(:,2),Morph(:,3),10,Morph(:,4));
%
%   Reference:
%   Betzel, Medaglia, Bassett (2018). Diversity of meso-scale architecture
%   in human and non-human connectomes. Nature Communications, 9(1), 346.
%
%   Richard Betzel, University of Pennsylvania, 2018
%

if nargin == 2
    verbose = false;
end

[N,NReps] = size(Ci);
C = zeros(N,5,NReps);
Morph = [];
for iRep = 1:NReps
    ci = Ci(:,iRep);
    I = dummyvar(ci);
    S = sum(I);
    
    if any(S == 1) & verbose
        warning('There are singleton communities - removing them');
    end
    
    D = (S'*S) - diag(S);
    
    AD = (I'*(A*I))./D;
    rmv = find(S == 1);
    rmvNodes = zeros(length(rmv),1);
    for iNode = 1:length(rmv)
        rmvNodes(iNode) = find(ci == rmv(iNode));
    end
    AD(rmv,:) = [];
    AD(:,rmv) = [];
    I(:,rmv) = [];
    [~,K] = size(I);
    for i = 1:K
        for j = 1:K
            if i ~= j
                lab = 0;
                ADij = AD([i,j],[i,j]);
                % if minimum within community density > between, then motif is
                % assortative
                if min(ADij(1,1),ADij(2,2)) > ADij(1,2)
                    C(ci == i,1,iRep) = C(ci == i,1,iRep) + 1;
                    C(ci == j,1,iRep) = C(ci == j,1,iRep) + 1;
                    lab = 1;
                    % if maximum within community density < between, then motif is
                    % dis-assortative
                elseif max(ADij(1,1),ADij(2,2)) < ADij(1,2)
                    C(ci == i,2,iRep) = C(ci == i,2,iRep) + 1;
                    C(ci == j,2,iRep) = C(ci == j,2,iRep) + 1;
                    lab = 2;
                    % if within i is greater than between i,j and between i,j is
                    % greater than within j, then i is core, j is periphery
                elseif ADij(1,1) > ADij(1,2) && ADij(1,2) > ADij(2,2)
                    C(ci == i,3,iRep) = C(ci == i,3,iRep) + 1;
                    C(ci == j,4,iRep) = C(ci == j,4,iRep) + 1;
                    lab = 3;
                    % if within j is greater than between i,j and between i,j is
                    % greater than within i, then j is core, i is periphery
                elseif ADij(2,2) > ADij(1,2) && ADij(1,2) > ADij(1,1)
                    C(ci == j,3,iRep) = C(ci == j,3,iRep) + 1;
                    C(ci == i,4,iRep) = C(ci == i,4,iRep) + 1;
                    lab = 3;
                end
                Morph = [Morph; ADij(1,1), ADij(2,2), ADij(1,2),lab,i,j];
            end
        end
    end
    C(rmvNodes,:,iRep) = nan;
end
Prob = bsxfun(@rdivide,C(:,1:4,:),sum(C(:,1:4,:),2));
C(:,1:4,:) = Prob;
C(:,5,:) = -nansum(Prob.*log2(Prob),2);