%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TVB parameter space exploration + checks on modeling parameters
% 
% Written by Hannelore Aerts 
% (UGent, Faculty of Psychology, Department of Data Analysis)
% Date last modification: 05/02/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Settings
scale = 68;
thr = {'thrA', 'thrR'};
dist_metric = {'Len', 'Log'};

%% Loop over all subjects
main='/home/hannelore/Documents/ANALYSES/TVB_global2/subjects/';
sublist=dir(main); 
n=length(sublist);

for sub=1:n
    for tt = 1:2
        tm = thr{tt};
        for dd = 1:2
            dm = dist_metric{dd};
            
            
            %% Load empirical data
            
            subID = sublist(sub).name;
            disp(['Processing ' subID ' (scale' num2str(scale) ', ' tm ', dist' dm ')']);
            
            in_fold = fullfile('/home/hannelore/Documents/ANALYSES/BTC_prepro/subjects/preop', subID);
            out_fold = (['/home/hannelore/Documents/ANALYSES/TVB_global2/subjects/' subID '/TVBii_' subID '_G-Ji']);

            cd(fullfile(in_fold, 'fmri'));
            load([subID '_fMRI_new.mat']);
            FC_emp = FC_cc_DK68;
            FC_emp = weight_conversion(FC_emp, 'autofix'); %put diagonal to zero
            FC_emp_z=atanh(FC_emp); %Fisher Z transformation
            clear CON* PAT* FC_cc* FC_mi ROI*

            cd(out_fold);
            SCtmp=load(['TVBiiInput_' tm '_scale' num2str(scale) '_dist' dm '.mat'], ['SC' tm 'n']);
            SC=cell2mat(struct2cell(SCtmp));
            clear SCtmp


            %% Compare empirical structure vs. function
            
            upperIdx = ones([scale scale]);
            upperIdx = triu(upperIdx,1);
            len=length(find(upperIdx));

            FC_emp_zv = reshape(FC_emp_z(~~upperIdx), [len 1]);
            SCv = reshape(SC(~~upperIdx), [len 1]);

            SC_FC_corP = corr(SCv, FC_emp_zv);
            SC_FC_corS = corr(SCv, FC_emp_zv, 'Type', 'Spearman');


            %% Load simulated BOLD 

            params_path=out_fold;
            sim_path=([params_path '/output_' tm '_scale' num2str(scale) '_dist' dm]);
            paramslist=dir([params_path '/params/param_set_*']);
            m=length(paramslist);

            TVBii_PSE = nan(m,10);
            TVBii_Ji = nan(scale,m);

            for index=2:m
                % Load parameters
                params = dlmread([params_path '/params/param_set_' num2str(index)]);

                % Load simulated TS
                s=dir([sim_path '/BOLD_param_set_' num2str(index) '.txt']);

                if exist([sim_path '/BOLD_param_set_' num2str(index) '.txt']) == 0
                    disp('No file found: putting nan in output')
                    SC_FC_corP(2)=nan;
                    emp_sim_cor(2)=nan;
                    Ji = nan(68,1);
                    J = mean(Ji);

                elseif s.bytes ~= 0
                    Ji = dlmread([sim_path '/BOLD_param_set_' num2str(index) '.txt']);
                    Ji = Ji(1:scale,2);
                    J = mean(Ji);

                    TS_sim = dlmread([sim_path '/BOLD_param_set_' num2str(index) '.txt']);
                    TS_sim = TS_sim(end-179:end,:);
                    FC_sim = corr(TS_sim);
                    FC_sim = weight_conversion(FC_sim, 'autofix');
                    FC_sim_z = atanh(FC_sim);
                    FC_sim_zv = reshape(FC_sim_z(~~upperIdx), [len 1]);

                    % Correlate simulated and empirical BOLD
                    emp_sim_corP = corr(FC_emp_zv, FC_sim_zv);
                    emp_sim_corS = corr(FC_emp_zv, FC_sim_zv, 'Type', 'Spearman');


                elseif s.bytes == 0
                    disp('Empty file found: putting nan in output')
                    emp_sim_cor=nan;
                    Ji = nan(68,1);
                    J = mean(Ji);

                end;

                % Save output
                TVBii_PSE(index, 1) = params(1); %number of nodes
                TVBii_PSE(index, 2) = params(2); %G
                TVBii_PSE(index, 3) = params(3); %J_NMDA
                TVBii_PSE(index, 4) = params(4); %W+
                TVBii_PSE(index, 5) = J; %average J
                TVBii_PSE(index, 6) = params(9); %TR
                TVBii_PSE(index, 7) = SC_FC_corP; %SC_FC correlation (Pearson)
                TVBii_PSE(index, 8) = SC_FC_corS; %SC_FC correlation (Spearman)
                TVBii_PSE(index, 9) = emp_sim_corP; %emp_sim cor (Pearson)
                TVBii_PSE(index, 10) = emp_sim_corS; %emp_sim cor (Spearman)

                TVBii_Ji(:,index)=Ji;
                
                clear FC_sim* TS_sim J Ji emp_sim_cor* params s
                
            end 
            
            cd(out_fold);
            save(sprintf('TVBiiOutput_%s_scale%d_dist%s',tm,scale,dm), 'FC_emp', 'FC_emp_z', 'FC_emp_zv', 'SC', 'SCv', 'TVBii_PSE', 'TVBii_Ji');
            clear TVBii_PSE TVBii_Ji params* sim_path SC_FC_cor* FC* SC* 
            clear in_fold index len out_fold upperIdx
            
        end
    end 
end 



%% Exploration of results: effect of thresholding, distance and correlation metric on G

results_fold='/home/hannelore/Documents/ANALYSES/TVB_global2/results_TVBii';

for sub=1:n
    subID = sublist(sub).name;
    disp(['Processing ' subID]);      
    out_fold = (['/home/hannelore/Documents/ANALYSES/TVB_global2/subjects/' subID '/TVBii_' subID '_G-Ji']);
    cd(out_fold);

    % Load all data
    load([out_fold '/TVBiiOutput_thrA_scale68_distLen'], 'TVBii_PSE')
    G_ALen = TVBii_PSE(:,2);
    corP_ALen = TVBii_PSE(:,9);
    corS_ALen = TVBii_PSE(:,10);
    clear TVBii_PSE
    
    load([out_fold '/TVBiiOutput_thrA_scale68_distLog'], 'TVBii_PSE')
    G_ALog = TVBii_PSE(:,2);
    corP_ALog = TVBii_PSE(:,9);
    corS_ALog = TVBii_PSE(:,10);
    clear TVBii_PSE
    
    load([out_fold '/TVBiiOutput_thrR_scale68_distLen'], 'TVBii_PSE')
    G_RLen = TVBii_PSE(:,2);
    corP_RLen = TVBii_PSE(:,9);
    corS_RLen = TVBii_PSE(:,10);
    clear TVBii_PSE
    
    load([out_fold '/TVBiiOutput_thrR_scale68_distLog'], 'TVBii_PSE')
    G_RLog = TVBii_PSE(:,2);
    corP_RLog = TVBii_PSE(:,9);
    corS_RLog = TVBii_PSE(:,10);
    clear TVBii_PSE
    
    % Plot   
    cd(results_fold);
    
    % G ~ corP for different thresholds and distance metrics
    plot(G_ALen, corP_ALen, 'LineWidth', 2, 'Color', 'y');
    xlabel('G'); ylabel('Pearson correlation FCsim - FCemp');
    hold on
    plot(G_ALog, corP_ALog, 'LineWidth', 2, 'Color', 'b'); 
        % complete overlap with corP_ALen
    hold on
    plot(G_RLen, corP_RLen, 'LineWidth', 2, 'Color', 'c');
    hold on
    plot(G_RLog, corP_RLog, 'LineWidth', 2, 'Color', 'g'); 
        % complete overlap with corP_distLen
    legend('thrA, distLog', 'thrA, distLen', 'thrR, distLog', 'thrR, distLen');
    print([subID '_G-corP_thr-dist'], '-dpng');
    close all
    
    % G ~ corS for different thresholds and distance metrics
    plot(G_ALen, corS_ALen, 'LineWidth', 2, 'Color', 'y');
    xlabel('G'); ylabel('Pearson correlation FCsim - FCemp');
    hold on
    plot(G_ALog, corS_ALog, 'LineWidth', 2, 'Color', 'b');
        % complete overlap with corS_ALen
    hold on
    plot(G_RLen, corS_RLen, 'LineWidth', 2, 'Color', 'c');
    hold on
    plot(G_RLog, corS_RLog, 'LineWidth', 2, 'Color', 'g'); 
        % complete overlap with corS_RLen
    legend('thrA, distLog', 'thrA, distLen', 'thrR, distLog', 'thrR, distLen');
    print([subID '_G-corS_thr-dist'], '-dpng');
    close all
    
    % G ~ cor for different thresholds and correlation metrics
    plot(G_ALog, corP_ALog, 'LineWidth', 2, 'Color', 'y');
    xlabel('G'); ylabel('Correlation FCsim - FCemp');
    hold on
    plot(G_ALog, corS_ALog, 'LineWidth', 2, 'Color', 'b');
    hold on
    plot(G_RLog, corP_RLog, 'LineWidth', 2, 'Color', 'c');
    hold on
    plot(G_RLog, corS_RLog, 'LineWidth', 2, 'Color', 'g'); 
    legend('thrA, corP', 'thrA, corS', 'thrR, corP', 'thrR, corS');
    print([subID '_G-cor_thr-cor'], '-dpng');
    close all
        % close overlap between threshold methods, greatest difference
        % between correlation metrics but qualitatively same pattern
        % (though Pearson cor generally higher than Spearman cor)
    
    clear G_* cor*
end

    %--> Continue with distLog, corP, both thresholds!

    
%% Save Gmax, SC-FC, FCemp-FCsim, and J for both thresholds

results_fold='/home/hannelore/Documents/ANALYSES/TVB_global2/results_TVBii';
TVBii_results_thrA=nan(n,5);
TVBii_results_thrR=nan(n,5);

for sub=1:n
    subID = sublist(sub).name;
    disp(['Processing ' subID]);      
    out_fold = (['/home/hannelore/Documents/ANALYSES/TVB_global2/subjects/' subID '/TVBii_' subID '_G-Ji']);
    cd(out_fold);

    % ThrA
    load([out_fold '/TVBiiOutput_thrA_scale68_distLog'], 'TVBii_PSE');
    results_sorted = sortrows(TVBii_PSE, -9);

    TVBii_results_thrA(sub, 1)=sub; %subject number (alphabetical order)
    TVBii_results_thrA(sub, 2)=results_sorted(1,7); %SC-FC
    TVBii_results_thrA(sub, 3)=results_sorted(1,9); %corP FCemp-FCsim max
    TVBii_results_thrA(sub, 4)=results_sorted(1,2); %G
    TVBii_results_thrA(sub, 5)=results_sorted(1,5); %average J  
    
    clear TVBii_PSE results_sorted
    
    % ThrR
    load([out_fold '/TVBiiOutput_thrR_scale68_distLog'], 'TVBii_PSE');
    results_sorted = sortrows(TVBii_PSE, -9);

    TVBii_results_thrR(sub, 1)=sub; %subject number (alphabetical order)
    TVBii_results_thrR(sub, 2)=results_sorted(1,7); %SC-FC
    TVBii_results_thrR(sub, 3)=results_sorted(1,9); %corP FCemp-FCsim max
    TVBii_results_thrR(sub, 4)=results_sorted(1,2); %G
    TVBii_results_thrR(sub, 5)=results_sorted(1,5); %average J  
    
    clear TVBii_PSE results_sorted
end 
    
% Save
cd(results_fold)
save('results_G-J_thrA', 'TVBii_results_thrA');
save('results_G-J_thrR', 'TVBii_results_thrR');




%% Sanity check 1: cor FCemp-FCsim using CONavg parameters

for tt = 1:2
    tm = thr{tt};
    dm = 'Log';
    for sub=1:n
        % Load empirical data
        subID = sublist(sub).name;
        disp(['Processing ' subID ' (scale' num2str(scale) ', ' tm ', dist' dm ')']);
            
        in_fold = fullfile('/home/hannelore/Documents/ANALYSES/BTC_prepro/subjects/preop', subID);
        out_fold = (['/home/hannelore/Documents/ANALYSES/TVB_global2/subjects/' subID '/TVBii_' subID '_G-Ji']);

        cd(fullfile(in_fold, 'fmri'));
        load([subID '_fMRI_new.mat']);
        FC_emp = FC_cc_DK68;
        FC_emp = weight_conversion(FC_emp, 'autofix'); %put diagonal to zero
        FC_emp_z=atanh(FC_emp); #Fisher-Z transform
        clear CON* PAT* FC_cc* FC_mi ROI*
        upperIdx = ones([scale scale]);
        upperIdx = triu(upperIdx,1);
        len=length(find(upperIdx));
        FC_emp_zv = reshape(FC_emp_z(~~upperIdx), [len 1]);

        % Load simulated BOLD 
        sim_path=([out_fold '/output_' tm '_scale' num2str(scale) '_dist' dm '_CONavg']);
        TS_sim = dlmread([sim_path '/BOLD_param_set_CON.txt']);
        TS_sim = TS_sim(end-179:end,:);
        FC_sim = corr(TS_sim);
        FC_sim = weight_conversion(FC_sim, 'autofix');
        FC_sim_z = atanh(FC_sim);
        FC_sim_zv = reshape(FC_sim_z(~~upperIdx), [len 1]);
        clear upperIdx len

        % Correlate simulated and empirical BOLD
        emp_sim_corP(sub,1) = corr(FC_emp_zv, FC_sim_zv);
        emp_sim_corS(sub,1) = corr(FC_emp_zv, FC_sim_zv, 'Type', 'Spearman');

        clear FC_* TS_sim in_fold out_fold
    end 
    results_fold='/home/hannelore/Documents/ANALYSES/TVB_global2/results_TVBii';
    cd(results_fold)
    emp_sim_corP_CONavgparams=emp_sim_corP;
    emp_sim_corS_CONavgparams=emp_sim_corS;
    save(['SanChecks_' tm], 'emp_sim_corP_CONavgparams', 'emp_sim_corS_CONavgparams');
end

% Plot subject max FCemp-FCsim correlation against that obtained with
% average parameters

% thrA
cd(results_fold)
load('SanChecks_thrA.mat');
load('results_G-J_thrA.mat');
dif = TVBii_results_thrA(:,3) - emp_sim_corP_CONavgparams;

% Difference between 2 distributions significant?
[h,p,ks2stat] = kstest2(emp_sim_corP_CONavgparams, TVBii_results_thrA(:,3))
    % --> yes! ks2stat=0.3889, p=0.0059
clear emp_sim_cor* TVBii_results_thrA dif h p ks2stat
          

% thrR
cd(results_fold)
load('SanChecks_thrR.mat');
load('results_G-J_thrR.mat');
dif = TVBii_results_thrR(:,3) - emp_sim_corP_CONavgparams;

% Difference between 2 distributions significant?
[h,p,ks2stat] = kstest2(emp_sim_corP_CONavgparams, TVBii_results_thrR(:,3))
    % --> yes! ks2stat=0.3889, p=0.0059
clear emp_sim_cor* TVBii_results_thrA dif h p ks2stat



%% Sanity check 2: use average SC matrix with subject-specific model params
% (optimized G + median J [because of skewnes])

% First get CON average SC weights & dist matrix
tm = 'thrA';
dm = 'Log';
for sub=1:11
    % Load individual SC matrix and add it to 3D array
    subID = sublist(sub).name;
    disp(['Processing ' subID ' (scale' num2str(scale) ', ' tm ', dist' dm ')']);
            
    out_fold = (['/home/hannelore/Documents/ANALYSES/TVB_global2/subjects/' subID '/TVBii_' subID '_G-Ji']);
    cd(out_fold)
    load(['TVBiiOutput_' tm '_scale' num2str(scale) '_dist' dm '.mat'], 'SC');
    SCall(:,:,sub)=SC;
    clear SC
    load(['TVBiiInput_' tm '_scale' num2str(scale) '_dist' dm '.mat'], 'SC_dist');
    SCdistall(:,:,sub)=SC_dist;
    clear SC_dist
end
SC_CONavg=mean(SCall, 3);
SCdist_CONavg=mean(SCdistall, 3);
    %now run "Generate_TVBii_Input_v4.m"!
CONavg_fold=([main 'avgCON']);
cd(CONavg_fold)
save('SC_CONavg', 'SC_CONavg', 'SCdist_CONavg');  


% Compare FCsim (using SCavg) to FCemp 
main='/home/hannelore/Documents/ANALYSES/TVB_global2/subjects/';
sublist=dir(main); sublist=sublist(3:38,:);
n=length(sublist);

for sub=1:n
    % Load empirical data
    subID = sublist(sub).name;
    disp(['Processing ' subID ' (scale' num2str(scale) ', ' tm ', dist' dm ')']);     
    in_fold = ([main subID '/TVBii_' subID '_G-Ji']);
    load([in_fold '/TVBiiOutput_thrA_scale68_distLen.mat'], 'FC_emp_zv'); 

    % Load simulated BOLD 
    sim_path=([main 'avgCON/TVBii_G-Ji/output']);
    TS_sim = dlmread([sim_path '/BOLD_param_set_' num2str(sub) '.txt']);
    TS_sim = TS_sim(end-179:end,:);
    FC_sim = corr(TS_sim);
    FC_sim = weight_conversion(FC_sim, 'autofix');
    FC_sim_z = atanh(FC_sim);
    upperIdx=ones([scale scale]);upperIdx=triu(upperIdx,1);len=length(find(upperIdx));
    FC_sim_zv = reshape(FC_sim_z(~~upperIdx), [len 1]);
    clear upperIdx len

    % Correlate simulated and empirical BOLD
    emp_sim_corP(sub,1) = corr(FC_emp_zv, FC_sim_zv);
    emp_sim_corS(sub,1) = corr(FC_emp_zv, FC_sim_zv, 'Type', 'Spearman');

    clear FC_* TS_sim in_fold out_fold
end

results_fold='/home/hannelore/Documents/ANALYSES/TVB_global2/results_TVBii';
cd(results_fold)
emp_sim_corP_CONavgSC=emp_sim_corP;
emp_sim_corS_CONavgSC=emp_sim_corS;
save('SanChecks_thrA', 'emp_sim_corP_CONavgSC', 'emp_sim_corS_CONavgSC', '-append');
    

% Difference between 2 distributions significant?
[h,p,ks2stat] = kstest2(emp_sim_corP_CONavgSC, TVBii_results_thrA(:,3))
    % --> yes! ks2stat=0.4444, p=0.001 

   

%% Sanity check 3: use average SC matrix to optimize indiv modelparams 

% First get CON average SC weights & dist matrix
cd([TVB_fold 'avgCON']);
load('SC_CONavg', 'SC_CONavg');  

% Compare FCsim (using SCavg) to FCemp 
main='/home/hannelore/Documents/ANALYSES/TVB_global2/subjects/';
sublist=dir(main); 
n=length(sublist);
upperIdx=ones([scale scale]);upperIdx=triu(upperIdx,1);len=length(find(upperIdx));

% Initialise output 
emp_sim_cor=nan(n,200);

for sub=1:n
    % Load empirical data
    subID = sublist(sub).name;
    disp(['Processing ' subID]);     
    in_fold = ([main subID '/TVBii_' subID '_G-Ji']);
    load([in_fold '/TVBiiOutput_thrA_scale68_distLen.mat'], 'FC_emp_zv'); 
    
    % Load simulated BOLD 
    sim_path=([main 'avgCON/TVBii_G-Ji/output2']);
    params_path=([main 'avgCON/TVBii_G-Ji/params2']);
    paramslist=dir([params_path '/param_set_*']);
    m=length(paramslist);
            
    for index=1:m
        % For now, only save correlation between FCsim & FCemp
        TS_sim = dlmread([sim_path '/BOLD_param_set_' num2str(index) '.txt']);
        TS_sim = TS_sim(end-179:end,:);
        FC_sim = corr(TS_sim);
        FC_sim = weight_conversion(FC_sim, 'autofix');
        FC_sim_z = atanh(FC_sim);
        FC_sim_zv = reshape(FC_sim_z(~~upperIdx), [len 1]);

        % Correlate simulated and empirical BOLD
        emp_sim_cor(sub, index) = corr(FC_emp_zv, FC_sim_zv);
        
        % Clean up
        clear TS* FC_sim*
    end
    clear FC_emp_zv
end

% Save output
output_fold=('/home/hannelore/Documents/ANALYSES/TVB_global2/results_TVBii');
cd(output_fold);
emp_sim_corP_CONavgSCparams=emp_sim_cor;
save('SanChecks_thrA.mat', 'emp_sim_corP_CONavgSCparams', '-append');
emp_sim_corP_CONavgSCparams_max = max(emp_sim_corP_CONavgSCparams');
save('SanChecks_thrA.mat', 'emp_sim_corP_CONavgSCparams_max', '-append');                


% Difference between 2 distributions significant?
load([output_fold '/results_G-J_thrA.mat']);
[h,p,ks2stat] = kstest2(emp_sim_corP_CONavgSCparams_max, TVBii_results_thrA(:,3))
    % --> not better, nor worse! ks2stat=0.1111, p=0.97


% Plot for couple of subjects distribution of G-cor to see whether same G
% values are obtained with both methods

sub='CON03T1' %n=3
sub='CON09T1' %n=9
sub='PAT03T1' %n=14
sub='PAT07T1' %n=18

load([main '/' sub '/TVBii_' sub '_G-Ji/TVBiiOutput_thrA_scale68_distLog.mat'], 'TVBii_PSE');
plot(TVBii_PSE(:,9), 'Linewidth', 2, 'Color', 'b')
hold on
plot(emp_sim_corP_CONavgSCparams(18,:), 'Linewidth', 2, 'Color', 'c')
legend('Subject-specific params & SC', 'CON avg SC + PSE') 
legend('Location', 'southoutside'); legend('boxoff')
xlabel('G'); ylabel('Pearson correlation coefficient');
cd([output_fold '/individual plots'])
print([sub '_indivSC-avgSCpse'], '-dpng')





%% Check Ji distribution : skewed

J_avg=nan(n,1);
J_md=nan(n,1);

for tt = 2
    tm = thr{tt};
        
    for index=1:n
        dm = 'Log';  
        scale = 68;

        % Load data
        subID=sublist(index).name;
        disp(['Processing ' subID ' (scale' num2str(scale) ', ' tm ', dist' dm ')']);
        load([main subID '/TVBii_' subID '_G-Ji/TVBiiOutput_' tm '_scale' num2str(scale) '_dist' dm '.mat'], 'Ji_max')

        cd('/home/hannelore/Documents/ANALYSES/TVB_global2/results_TVBii/individual plots')
        figure();hist(Ji_max)
        %print([subID '_Ji_dist'], '-dpng');
        
        J_avg(index,1)=mean(Ji_max);
        J_md(index,1)=median(Ji_max);
        
        clear Ji_max
    end
end



%% Average Ji in vicinity of tumor

% Tumor nodes
DM_fold='/home/hannelore/Documents/ANALYSES/BTC_prepro/tumorrois_all.csv';
TVB_fold='/home/hannelore/Documents/ANALYSES/TVB_global2/subjects/';
tumornodes=csvread(DM_fold, 1);

cd(TVB_fold);
sublist=dir(TVB_fold); sublist=sublist(3:end,:);
n=length(sublist);

for tt = 1:2
    tm = thr{tt};
    for sub=2:n
        dm = 'Log';
        
        % First get Ji values corresponding to max FCsim-FCemp
        subID = sublist(sub).name;
        disp(['Processing ' subID ' (scale' num2str(scale) ', ' tm ', dist' dm ')']);
        in_path = [TVB_fold subID '/TVBii_' subID '_G-Ji/'];
        cd(in_path);
        load([in_path 'TVBiiOutput_' tm '_scale' num2str(scale) '_dist' dm '.mat'], 'TVBii_PSE', 'TVBii_Ji');
        
        [maxval, maxind] = max(TVBii_PSE(:,9));
        Ji_max = TVBii_Ji(:,maxind);
        
        % Now identify tumor nodes
        sub_tumornodes = find(tumornodes(:,sub));
        
        Ji_tumornodes = Ji_max(sub_tumornodes,:);
        J_tumornodes = mean(Ji_tumornodes);
        J_tumornodes_allsubs(sub,1) = mean(Ji_tumornodes); 
        J_tumornodes_allsubs_md(sub,1) = median(Ji_tumornodes); 
        %save([in_path 'TVBiiOutput_' tm '_scale' num2str(scale) '_dist' dm '.mat'], 'J_tumornodes', 'Ji_tumornodes', 'Ji_max', '-append'); 
        clear subID in_path Ji_* J_tumornodes sub_tumornodes TVBii_PSE TVBii_Ji maxval maxind
    end
    clear J_tumornodes_allsubs*
end
clear tumornodes DM_fold


% Tumor nodes (random nodes for CON instead of MEN nodes)

DM_fold='/home/hannelore/Documents/ANALYSES/BTC_prepro/tumorrois_all_CONrandom.csv';
TVB_fold='/home/hannelore/Documents/ANALYSES/TVB_global2/subjects/';
tumornodes=csvread(DM_fold, 1);
numberofrois=sum(tumornodes~=0,1)'

cd(TVB_fold);
sublist=dir(TVB_fold); sublist=sublist(3:end,:);
n=length(sublist);

for tt = 1:2
    tm = thr{tt};
    for sub=1:n
        dm = 'Log';
        
        % First get Ji values corresponding to max FCsim-FCemp
        subID = sublist(sub).name;
        disp(['Processing ' subID ' (scale' num2str(scale) ', ' tm ', dist' dm ')']);
        in_path = [TVB_fold subID '/TVBii_' subID '_G-Ji/'];
        cd(in_path);
        load([in_path 'TVBiiOutput_' tm '_scale' num2str(scale) '_dist' dm '.mat'], 'TVBii_PSE', 'TVBii_Ji');
        
        [maxval, maxind] = max(TVBii_PSE(:,9));
        Ji_max = TVBii_Ji(:,maxind);
        
        % Now identify tumor nodes
        sub_tumornodes = find(tumornodes(:,sub));
        
        Ji_tumornodes = Ji_max(sub_tumornodes,:);
        J_tumornodes = mean(Ji_tumornodes);
        J_tumornodes_allsubs(sub,1) = mean(Ji_tumornodes); 
        J_tumornodes_allsubs_md(sub,1) = median(Ji_tumornodes); 
        %save([in_path 'TVBiiOutput_' tm '_scale' num2str(scale) '_dist' dm '.mat'], 'J_tumornodes', 'Ji_tumornodes', 'Ji_max', '-append'); 
        clear subID in_path Ji_* J_tumornodes sub_tumornodes TVBii_PSE TVBii_Ji maxval maxind
    end
    clear J_tumornodes_allsubs*
end
clear tumornodes DM_fold



% Check median Ji in non-tumor regions

DM_fold='/home/hannelore/Documents/ANALYSES/BTC_prepro/tumorrois_all_CONrandom.csv';
TVB_fold='/home/hannelore/Documents/ANALYSES/TVB_global2/subjects/';
tumornodes=csvread(DM_fold, 1);
numberofrois=sum(tumornodes~=0,1)'

cd(TVB_fold);
sublist=dir(TVB_fold); sublist=sublist(3:end,:);
n=length(sublist);

for tt = 1:2
    tm = thr{tt};
    for sub=1:n
        dm = 'Log';
        
        % First get Ji values corresponding to max FCsim-FCemp
        subID = sublist(sub).name;
        disp(['Processing ' subID ' (scale' num2str(scale) ', ' tm ', dist' dm ')']);
        in_path = [TVB_fold subID '/TVBii_' subID '_G-Ji/'];
        cd(in_path);
        load([in_path 'TVBiiOutput_' tm '_scale' num2str(scale) '_dist' dm '.mat'], 'TVBii_PSE', 'TVBii_Ji');
        
        [maxval, maxind] = max(TVBii_PSE(:,9));
        Ji_max = TVBii_Ji(:,maxind);
        
        % Now identify tumor nodes
        sub_nontumornodes = find(tumornodes(:,sub)==0);
        Ji_nontumornodes = Ji_max(sub_nontumornodes,:);
        J_nontumornodes_allsubs_md(sub,1) = median(Ji_nontumornodes); 
        clear subID in_path Ji_* sub_nontumornodes TVBii_PSE TVBii_Ji maxval maxind
    end
    clear J_tumornodes_allsubs*
end
clear tumornodes DM_fold



%% Check group differences in Ji_tumor persist when regressing out SC strength

main='/home/hannelore/Documents/ANALYSES/TVB_global2/subjects/';
sublist=dir(main); sublist=sublist(3:end,:);
n=length(sublist);
DM_fold='/home/hannelore/Documents/ANALYSES/BTC_prepro/tumorrois_all_CONrandom.csv';
tumornodes=csvread(DM_fold, 1);
tm='thrA'
dm='Log'

for sub=1:n
   subID = sublist(sub).name;
   disp(['Processing ' subID ' (scale' num2str(scale) ', ' tm ', dist' dm ')']);
   
   % Load Ji data corresponding to optimal G value
   load([main subID '/TVBii_' subID '_G-Ji/TVBiiOutput_thrA_scale68_distLog.mat'], 'TVBii_PSE');
   [maxcor, maxcorInd] = max(TVBii_PSE(:,9));

   % Load Ji and SC strength values for corresponding max cor value of G
   results_fold=([main subID '/TVBii_' subID '_G-Ji/output_thrA_scale68_distLog']);
   J = dlmread([results_fold '/BOLD_param_set_' num2str(maxcorInd) '.txt']); 
   Ji = J(1:scale,2);
   Ji_strength = J(1:scale,1);
   B=regress(Ji, Ji_strength);
   res(:,sub)=Ji-B*Ji_strength;   
   
   % Median JiBrain res
   J_res(sub)=median(res(:,sub));
   
   % Median JiTumor res
   sub_tumornodes = find(tumornodes(:,sub));
   JiT_res(sub) = median(res(sub_tumornodes,sub));
   
   % Median JiNonTumor res
   sub_nontumornodes = find(tumornodes(:,sub)==0);
   JiNT_res(sub) = median(res(sub_nontumornodes,sub));
   
   clear J TVBii_PSE B Ji Ji_strength
end



%% Check group differences in Ji_tumor persist when regressing out only ROI size

for sub=1:n
   subID = sublist(sub).name;
   disp(['Processing ' subID ' (scale' num2str(scale) ', ' tm ', dist' dm ')']);
   
   % Load Ji data corresponding to optimal G value
   load([main subID '/TVBii_' subID '_G-Ji/TVBiiOutput_thrA_scale68_distLog.mat'], 'TVBii_PSE');
   [maxcor, maxcorInd] = max(TVBii_PSE(:,9));

   % Load Ji and SC strength values for corresponding max cor value of G
   results_fold=([main subID '/TVBii_' subID '_G-Ji/output_thrA_scale68_distLog']);
   J = dlmread([results_fold '/BOLD_param_set_' num2str(maxcorInd) '.txt']); 
   Ji = J(1:scale,2);
   Ji_strength = J(1:scale,1);
   sub_roisize=RoiSize(:,sub);
   X = sub_roisize;
   
   [B,bint,res]=regress(Ji, X);
   
   % Median JiBrain res
   J_res3(sub)=median(res(:));
   
   % Median JiTumor res
   sub_tumornodes = find(tumornodes(:,sub));
   JiT_res3(sub) = median(res(sub_tumornodes,:));
   
   % Median JiNonTumor res
   sub_nontumornodes = find(tumornodes(:,sub)==0);
   JiNT_res3(sub) = median(res(sub_nontumornodes,:));
   
   clear maxcor* res results_fold J TVBii_PSE B Ji Ji_strength 
end


