%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GTA analyses using BCT on SC matrices.
% 
% Written by Hannelore Aerts 
% (UGent, Faculty of Psychology, Department of Data Analysis)
% Date last modification: 28/11/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_path = '/home/hannelore/Documents/ANALYSES/TVB_global2/subjects/';
cd(data_path);
sublist=dir('*');
m=length(sublist);

thr = {'thrA', 'thrR'};
for tt = 1:2
    tm = thr{tt};
        
    % Prepare output files
    GTA_results = nan(m,11);
    results_path = '/home/hannelore/Documents/ANALYSES/TVB_global2/results_GTA/';
    cd(results_path)
    fileID = fopen(sprintf('GTA_%s_results.txt', tm), 'w');
    fprintf(fileID, '%20s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s\n', 'subID', ['SC_' tm '_density'], ['SC_' tm '_clust'], ['SC_' tm '_Eloc'], ['SC_' tm '_Q'], ['SC_' tm '_Eglob'], ['SC_' tm '_comm'], ['SC_' tm '_degree'], ['SC_' tm '_strength'], ['SC_' tm '_EBC'], ['SC_' tm '_BC'], ['SC_' tm '_PC']); 
    cd(data_path)

    for index=1:m
        dm = 'Log';  
        scale = 68;

        %% Load data
        subID=sublist(index).name;
        disp(['Processing ' subID ' (scale' num2str(scale) ', ' tm ', dist' dm ')']);
        load([subID '/TVBii_' subID '_G-Ji/TVBiiOutput_' tm '_scale' num2str(scale) '_dist' dm '.mat'], 'SC')
        SC_dist=-log(SC); SC_dist(isinf(SC_dist)) = 0;

        %% Graph measures -- general
        % Density 
        density = density_und(SC);


        %% Graph measures -- segregation
        % Clustering coefficient (all weights must be between 0 and 1)
        clust = clustering_coef_wu(SC);
        clust_mean = mean(clust);

        % Local efficiency (all weights must be between 0 and 1)
        Eloc = efficiency_wei(SC, 1);
        Eloc_mean = mean(Eloc);

        % Modularity 
        M = nan(scale, 100);
        Q = nan(100,1);
        for i=1:100
            [M(:,i), Q(i,:)] = community_louvain(SC, [],[], 'modularity');
        end
        [Q_max, Q_ind] = max(Q);
        M_max = M(:,Q_ind);


        %% Graph measures -- integration
        % Global efficiency
        Eglob = efficiency_wei(SC);

        % Communicability
        comm = expm(SC); 
        comm_mean = mean(comm(:));


        %% Graph measures -- centrality
        % Degree
        degree = degrees_und(SC);
        degree_mean = mean(degree);

        % Strength
        strength = strengths_und(SC);
        strength_mean = mean(strength);

        % Betweenness centrality --> input is connection lengths instead
        % of connection weights
        [EBC, BC] = edge_betweenness_wei(SC_dist);
        EBC_mean = mean(EBC(:));
        BC_mean = mean(BC);

        % Participation coefficient
        PC = participation_coef(SC, M_max);
        PC_mean = mean(PC);


        %% Save results
        GTA_results(index, 1) = density;
        GTA_results(index, 2) = clust_mean;
        GTA_results(index, 3) = Eloc_mean;
        GTA_results(index, 4) = Q_max;
        GTA_results(index, 5) = Eglob;
        GTA_results(index, 6) = comm_mean;
        GTA_results(index, 7) = degree_mean;
        GTA_results(index, 8) = strength_mean;
        GTA_results(index, 9) = EBC_mean;
        GTA_results(index, 10) = BC_mean;
        GTA_results(index, 11) = PC_mean;

        fprintf(fileID, '%20s %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f\n', subID, GTA_results(index,:));

        clear BC* EBC* Eglob Eloc* M* PC* Q* clust* comm* degree* 
        clear density i strength* SC*    
        cd(data_path)
    end  
    cd(results_path);
    fclose(fileID);

    clear fileID index subID YeoOrder ans
    save(['GTA_' tm '_results_all']);
    clear GTA_results
end
    

