clc;clear;
cd E:\Github\hippocampus_cortex_gradient_youth % replace this with absolute path of your working directory
addpath(genpath('./Dependencies/Matlab/'))

% root path to save all the results
if ~exist('./Results','dir')
   mkdir('./Results') 
end
% path to save the results of hippocampal-cortical FCs and gradients
conn_indiv='./Results/conn_mat/';
avg_result='./Results/gradient_avg';
if ~exist(avg_result,'dir')
   mkdir(avg_result) 
end
if ~exist(conn_indiv,'dir')
   mkdir(conn_indiv) 
end

% path of the fMRI data 
hipp_fmri_dir='./Data/fMRI_data/hipp_native';
cortex_fmri_dir='./Data/fMRI_data/cortex_glasser';

% read the demographic information of HCP-D data
opts=detectImportOptions('./Code/HCD_LS_2.0_subject_completeness_bh.csv');
demographic_info=readtable('./Code/HCD_LS_2.0_subject_completeness_bh.csv',opts);
subjlist=demographic_info{:,1};
subj_num=length(subjlist);

hemis={'L','R'};

%initialize a cell to store the mean FC for left and right hippocampus 
conn_matrix_mean_2=cell(1,2); 
%initialize two cells to store the relevant results of gradient for left and right hippocampus 
gm_2=cell(1,2); 
gradient_2=cell(1,2);

%vertex number of hippocampal mid-thickness surface
vert_num=483;

if ~exist([avg_result '/avg-gradient-4-LR-' num2str(subj_num) '-dm.mat'],'file')

for h =1:2    
hemi=hemis{h};

conn_matrix=zeros(vert_num,360);
conn_matrix_z=zeros(vert_num,360);
k=0;

if isempty(dir([avg_result '/conn_matrix_mean_' hemi '-*.mat'])) 
parfor i=1:subj_num
subj_id=[subjlist{i} '_V1_MR'];
%disp(subj_id);

if ~exist([conn_indiv '/' subj_id '_' hemi '_conn_mat.mat'],'file')

if exist([hipp_fmri_dir '/' subj_id '/rfMRI_REST_hipp_ribbon_bandpass_native_hemi-' hemi '.func.gii'],'file') && ...
   exist([hipp_fmri_dir '/' subj_id '/rfMRI_REST_dentate_bandpass_native_hemi-' hemi '.func.gii'],'file') 

    fmri_cortex=cifti_read([cortex_fmri_dir '/' subj_id '/MNINonLinear/Results/rfMRI_REST/rfMRI_REST_Atlas_bp_360.ptseries.nii']);
    
    fmri_hipp=gifti([hipp_fmri_dir '/' subj_id '/rfMRI_REST_hipp_ribbon_bandpass_native_hemi-' hemi '.func.gii']);
    fmri_DG=gifti([hipp_fmri_dir '/' subj_id '/rfMRI_REST_dentate_bandpass_native_hemi-' hemi '.func.gii']);
    fmri_hippo=[fmri_hipp.cdata;fmri_DG.cdata];
else
    disp([subj_id,'- lacking fmri(hippo)'])
    continue;
end

C_r=corr(fmri_hippo',fmri_cortex.cdata');

if sum(sum(isnan(C_r))) ~=0
    disp(['the correlation matrix of ',subj_id,' has NAN']);
    if ~isempty(find(all(fmri_hippo==0,2), 1))
        disp(['the fmri data on the hippo-skeleton of ',subj_id,'-',hemi, ' has all 0 rows:'])
        disp(num2str(find(all(fmri_hippo==0,2))))
    end    
    if ~isempty(find(all(fmri_cortex.cdata==0,2), 1))
        disp(['the fmri data on the cerebral-cortex of ',subj_id,'-',hemi, ' has all 0 rows'])
    end   
    continue;
else
    % save the individual FCs
    parsave([conn_indiv '/' subj_id '_' hemi '_conn_mat.mat'], C_r)
    csvwrite([conn_indiv '/' subj_id '_' hemi '_conn_mat.csv'], C_r);
    C_z=atanh(C_r);
    conn_matrix_z=conn_matrix_z+C_z;
    k=k+1;

end

else
    % load individual FCs
    C_r=importdata([conn_indiv '/' subj_id '_' hemi '_conn_mat.mat']);
    C_z=atanh(C_r);
    conn_matrix_z=conn_matrix_z+C_z;
    k=k+1;  
end

end

% compute and save the mean FCs
conn_matrix_z_mean=conn_matrix_z/k;
conn_matrix_mean=tanh(conn_matrix_z_mean);
save([avg_result '/conn_matrix_mean_' hemi '-' num2str(k) '.mat'], 'conn_matrix_mean')
csvwrite([avg_result '/conn_matrix_mean_' hemi '-' num2str(k) '.csv'],conn_matrix_mean);

else
   % load mean FCs
   conn_struct=dir([avg_result '/conn_matrix_mean_' hemi '-*.mat']);
   load(strcat(avg_result, '/',conn_struct.name)); 
end
conn_matrix_mean_2{1}=conn_matrix_mean;

% Compute connectome gradient
if h==1
   % Build gradients using diffusion maps and normalized angle    
   gm = GradientMaps('kernel','na','approach','dm','n_components',10);
   % and fit to the data   
   gm = gm.fit(conn_matrix_mean);
   gradient_avg=gm.gradients{1,1};
else
     load([avg_result '/avg-gradient-4-L-' num2str(subj_num) '-dm.mat'])
    % Build gradients using diffusion maps and normalized angle
    % Align the gradients of right hippocampus to that of the left hippocampus
    gm = GradientMaps('kernel','na','approach','dm','n_components',10, 'alignment','pa'); 
    % and fit to the data
    gm = gm.fit(conn_matrix_mean, 'reference', gradient_avg);
    gradient_avg=gm.aligned{1,1};
end    

%save the results of group level hippocampal-cortical connectome gradients
%for each hemisphere
save([avg_result '/avg-gm-4-' hemi '-' num2str(subj_num) '-dm.mat'], 'gm')
gm_2{h}=gm;
save([avg_result '/avg-gradient-4-' hemi '-' num2str(subj_num) '-dm.mat'], 'gradient_avg')
gradient_2{h}=gradient_avg;
end

%save the results of group level hippocampal-cortical connectome gradients
%for both hemispheres
save([avg_result '/avg-gradient-4-LR-' num2str(subj_num) '-dm.mat'],'gradient_2');
save([avg_result '/avg-gm-4-LR-' num2str(subj_num) '-dm.mat'],'gm_2');

else

gradient_2=importdata([avg_result '/avg-gradient-4-LR-' num2str(subj_num) '-dm.mat']);
load([avg_result '/avg-gm-4-LR-' num2str(subj_num) '-dm.mat']);

end

%% Draw figures of group level gradients
% Path to save the surface related visualization
G_surf_figure=[avg_result '/G_surf_figure_avg_all'];
if ~exist(G_surf_figure,'dir')
   mkdir(G_surf_figure) 
end

% Path to save the scatter related visualization
G_scatter_figure=[avg_result '/G_scatter_figure_avg_all'];
if ~exist(G_scatter_figure,'dir')
   mkdir(G_scatter_figure) 
end

surf_temp='./Surf_temp'; 
spaces={'canonical', 'unfold'};

min_all=min([-gradient_2{1};-gradient_2{2}]);
max_all=max([-gradient_2{1};-gradient_2{2}]);

for h=1:2
hemi=hemis{h};
gradient_avg=-gradient_2{h};
gm=gm_2{h};

for sp=1:2
    space=spaces{sp};
    surf=SurfStatReadSurf([surf_temp '/tpl-avg_space-' space '_den-2mm_midthickness.obj']);
    if sp==1 && h==1
       surf.coord(1,:)=max(surf.coord(1,:))+0.5-surf.coord(1,:);
    elseif sp==2 && h==2
       surf.coord(2,:)=max(surf.coord(2,:))+0.5-surf.coord(2,:);
    end
    
   %% [Figure 1a] Illustrate the two mid-thickness surfaces (hipp, DG) with yellow and red corlor
    if h==1 && sp==1
        figure('Color','White','Units','Normalized','Position',[.1 .1 .56 .66]);
        [a,~]=BladeSurfStatViewData([ones(419,1);2*ones(64,1)], surf, ['Hipp-surf-', hemi, '_', space]);
        colormap([1,1,0;1,0,0])
        colorbar off
        saveas(a, [G_surf_figure '/Hipp-surf-', hemi, '_', space],'tif')       
    end    
    %% [Figure 1b] and [Extended Data Fig. 1]
    for g=1:5
        figure('Color','White','Units','Normalized','Position',[.1 .1 .56 .66]); 
        [a,~]=BladeSurfStatViewData(gradient_avg(:,g)', surf, ['Gradient1-', hemi, '_', space, '-avg']);
        colormap(viridis);camlight;SurfStatColLim([min_all(g) max_all(g)]);
        colorbar('ytick',linspace(min_all(g), max_all(g), 6))
        if sp==2;view(-90,90);end
        saveas(a, [G_surf_figure '/Gradient' num2str(g) '_',hemi, '_', space],'tif')
%         colorbar off
%         saveas(a, [G_surf_figure '/Gradient' num2str(g) '_',hemi, '_', space],'pdf')
    end    
end
    %% [Figure 1a] Draw figure of variance explained by each gradient component 
     f=scree_plot(gm.lambda{1});
     ExpVar=double(gm.lambda{1}./ sum(gm.lambda{1}));
     saveas(f.figure, [G_surf_figure '/ExpVar_',hemi],'pdf')
     
     text(double(1), ExpVar(1), ['  \leftarrow' num2str(ExpVar(1))]);
     text(double(2), ExpVar(2), ['  \leftarrow' num2str(ExpVar(2))]);
     text(double(3), ExpVar(3), ['  \leftarrow' num2str(ExpVar(3))]);
     text(double(4), ExpVar(4), ['  \leftarrow' num2str(ExpVar(4))]);
     text(double(5), ExpVar(5), ['  \leftarrow' num2str(ExpVar(5))]);
     saveas(f.figure, [G_surf_figure '/ExpVar_',hemi],'tif')
end

%% Test and visualize the relationships between hippocampal gradients and AP/PD coordinate
% load the pairwise vertex geodesic distance matrix of hippocampal mid-thickness surfaces
geo_dist_hipp=gifti([surf_temp '/hipp/midthickness_geoDist.gii']);geo_dist_hipp=geo_dist_hipp.cdata;
geo_dist_DG=gifti([surf_temp '/dentate/midthickness_geoDist.gii']);geo_dist_DG=geo_dist_DG.cdata;

%set the initial parameters for variogram algorithm
rand_init=0;n_surro=1000;
num_samp1=400;
num_neigh1=80;
num_samp2=60;
num_neigh2=10;

obj_subsample_hipp = variogram(geo_dist_hipp, 'ns', num_samp1, ...
    'knn', num_neigh1, 'random_state', rand_init, 'num_workers', 48);
obj_subsample_DG = variogram(geo_dist_DG, 'ns', num_samp2, ...
    'knn', num_neigh2, 'random_state', rand_init, 'num_workers', 48);

%% [Figure 1c] Test the relationships: G1-PA; G2-PA; G3-PD 
% store the p-values and r square values
p_all=nan(3,2);
r2_all=nan(3,2);

for h=1:2
    hemi=hemis{h};
    gradient_avg=-gradient_2{h};
    surf=SurfStatReadSurf([surf_temp '/tpl-avg_' hemi '_space-unfold_den-2mm_midthickness.obj']);
for g=1:3   
    if g==3
        if exist([surf_temp '/surrogates_subsample_PD_' hemi '.csv'], 'file')
            surrogates_subsample=csvread([surf_temp '/surrogates_subsample_PD_' hemi '.csv']);
        else
            surrogates_subsample_hipp = obj_subsample_hipp.fit(zscore(surf.coord(2,1:419))', n_surro);
            surrogates_subsample_DG = obj_subsample_DG.fit(zscore(surf.coord(2,420:483))', n_surro);
            surrogates_subsample=[surrogates_subsample_hipp;surrogates_subsample_DG];
            csvwrite([surf_temp '/surrogates_subsample_PD_' hemi '.csv'],surrogates_subsample)
        end
    else
        if exist([surf_temp '/surrogates_subsample_AP_' hemi '.csv'], 'file')
            surrogates_subsample=csvread([surf_temp '/surrogates_subsample_AP_' hemi '.csv']);
        else        
            surrogates_subsample_hipp = obj_subsample_hipp.fit(zscore(surf.coord(1,1:419))', n_surro);
            surrogates_subsample_DG = obj_subsample_DG.fit(zscore(surf.coord(1,420:483))', n_surro);
            surrogates_subsample=[surrogates_subsample_hipp;surrogates_subsample_DG];
            csvwrite([surf_temp '/surrogates_subsample_AP_' hemi '.csv'],surrogates_subsample)
        end
    end   
    
    if g==1
        r_real = corr(zscore(surf.coord(1,:))',gradient_avg(:,g));
        r_surrogate = corr(surrogates_subsample, gradient_avg(:,g));
        prctile_rank = mean(r_real > r_surrogate);
        r2_all(g,h)=r_real.^2;
    elseif g==2
        tbl=table(gradient_avg(:,g),zscore(surf.coord(1,:))','VariableNames',{'G2','AP'});
        lm=fitlm(tbl,'G2~1+AP+AP^2');
        Rsq_real = lm.Rsquared.Adjusted;
        r2_all(g,h)=Rsq_real;
        Rsq_surrogate=nan(n_surro,1);
        parfor i=1:n_surro            
            tbl=table(gradient_avg(:,g),surrogates_subsample(:,i),'VariableNames',{'G2','AP'});
            lm=fitlm(tbl,'G2~1+AP+AP^2');
            Rsq_surrogate(i) = lm.Rsquared.Adjusted;        
        end
        prctile_rank = mean(Rsq_real > Rsq_surrogate);
    else
        tbl=table(gradient_avg(:,g),zscore(surf.coord(2,:))','VariableNames',{'G3','PD'});
        lm=fitlm(tbl,'G3~1+PD+PD^2+PD^3');
        Rsq_real = lm.Rsquared.Adjusted;
        r2_all(g,h)=Rsq_real;
        Rsq_surrogate=nan(n_surro,1);
        parfor i=1:n_surro
            tbl=table(gradient_avg(:,g),surrogates_subsample(:,i),'VariableNames',{'G3','PD'});
            lm=fitlm(tbl,'G3~1+PD+PD^2+PD^3');
            Rsq_surrogate(i) = lm.Rsquared.Adjusted;        
        end
        prctile_rank = mean(Rsq_real > Rsq_surrogate);       
    end    
    
   if prctile_rank>0.5
      p_all(g,h)=1-prctile_rank;
   else
      p_all(g,h)=prctile_rank;
   end
   
%    if g==1
%        f=figure;h1=histfit(r_surrogate,100);set(h1(1),'facecolor','w');
%        hold on;plot([r_real,r_real],[0,20],'y-','LineWidth',4)
%        title(['G1(' hemi ')- AP coord: r=' num2str(r_real) ', p=' num2str(p_all(g,h))])
%        saveas(f, [G_scatter_figure '/G' num2str(g) '(' hemi ')_AP_coord'],'tif')    
%    elseif g==2
%        f=figure;h1=histfit(Rsq_surrogate,100);set(h1(1),'facecolor','w');
%        hold on;plot([Rsq_real,Rsq_real],[0,20],'y-','LineWidth',4)         
%        title(['G2(' hemi ')- AP coord: R-sq-adj=' num2str(Rsq_real) ', p=' num2str(p_all(g,h))])
%        saveas(f, [G_scatter_figure '/G' num2str(g) '(' hemi ')_AP_coord'],'tif')
%    else 
%        f=figure;h1=histfit(Rsq_surrogate,100);set(h1(1),'facecolor','w');
%        hold on;plot([Rsq_real,Rsq_real],[0,20],'y-','LineWidth',4)       
%        title(['G3(' hemi ')- PD coord: R-sq-adj=' num2str(Rsq_real) ', p=' num2str(p_all(g,h))])
%        saveas(f, [G_scatter_figure '/G' num2str(g) '(' hemi ')_PD_coord'],'tif')      
%    end
   
end
end

clear f
figure('Position',[100,100,1000,200]);
for h=1:2
    hemi=hemis{h};
    surf=SurfStatReadSurf([surf_temp '/tpl-avg_' hemi '_space-unfold_den-2mm_midthickness.obj']);   
    gradient_avg=-gradient_2{h};
for g=1:3 

    if g==1
        f(1,2*g+h-2)=gramm('x',zscore(surf.coord(1,:)), 'y',gradient_avg(:,g));
        f(1,2*g+h-2).set_names('x','z-scored P-A coordinate','y',['G' num2str(g) '-' hemi '-grp'])
        f(1,2*g+h-2).geom_point();f(1,2*g+h-2).stat_glm(); 
        f(1,2*g+h-2).set_point_options('markers',{'o'});
        if p_all(g,h)==0 
          f(1,2*g+h-2).set_title(['Rsq=' num2str(r2_all(g,h)) ', p<0.001']) 
        else
          f(1,2*g+h-2).set_title(['Rsq=' num2str(r2_all(g,h)) ', p=' num2str(p_all(g,h))])
        end
    elseif g==2
        f(1,2*g+h-2)=gramm('x',zscore(surf.coord(1,:)), 'y',gradient_avg(:,g));
        f(1,2*g+h-2).set_names('x','z-scored P-A coordinate','y',['G' num2str(g) '-' hemi '-grp'])
        f(1,2*g+h-2).geom_point();f(1,2*g+h-2).stat_fit('fun',@(b,c,d,x)b*x.^2+c*x+d,'intopt','functional');
        f(1,2*g+h-2).set_point_options('markers',{'o'});
        f(1,2*g+h-2).set_title(['Rsq=' num2str(r2_all(g,h)) ', p=' num2str(p_all(g,h))])
        if p_all(g,h)==0 
          f(1,2*g+h-2).set_title(['Rsq=' num2str(r2_all(g,h)) ', p<0.001']) 
        else
          f(1,2*g+h-2).set_title(['Rsq=' num2str(r2_all(g,h)) ', p=' num2str(p_all(g,h))])
        end
    else
        f(1,2*g+h-2)=gramm('x',zscore(surf.coord(2,:)), 'y',gradient_avg(:,g));
        f(1,2*g+h-2).set_names('x','z-scored P-D coordinate','y',['G' num2str(g) '-' hemi '-grp'])
        f(1,2*g+h-2).geom_point();f(1,2*g+h-2).stat_fit('fun',@(a,b,c,d,x)a*x.^3+b*x.^2+c*x+d,'intopt','functional');
        f(1,2*g+h-2).set_point_options('markers',{'o'});
        if p_all(g,h)==0 
          f(1,2*g+h-2).set_title(['Rsq=' num2str(r2_all(g,h)) ', p<0.001']) 
        else
          f(1,2*g+h-2).set_title(['Rsq=' num2str(r2_all(g,h)) ', p=' num2str(p_all(g,h))])
        end
    end      
    
end
end
f.set_title(['Visualization of the relationship between Gradient and AP(or PD) coordinate'])
f.set_color_options('map','brewer2')
f.draw();
f.export('file_name','G-corr-AP_PD','export_path',[G_scatter_figure '/'],'file_type','pdf','width',60,'height',10)

%% [Extended Data Fig. 2] Alternative relationships: G1-PD; G2-PD; G3-PA
p_all=nan(3,2);
r2_all=nan(3,2);
for h=1:2
    hemi=hemis{h};
    gradient_avg=-gradient_2{h};
    surf=SurfStatReadSurf([surf_temp '/tpl-avg_' hemi '_space-unfold_den-2mm_midthickness.obj']);
for g=1:3   
    if g==3
        if exist([surf_temp '/surrogates_subsample_AP_' hemi '.csv'], 'file')
            surrogates_subsample=csvread([surf_temp '/surrogates_subsample_AP_' hemi '.csv']);
        else        
            surrogates_subsample_hipp = obj_subsample_hipp.fit(zscore(surf.coord(1,1:419))', n_surro);
            surrogates_subsample_DG = obj_subsample_DG.fit(zscore(surf.coord(1,420:483))', n_surro);
            surrogates_subsample=[surrogates_subsample_hipp;surrogates_subsample_DG];
            csvwrite([surf_temp '/surrogates_subsample_AP_' hemi '.csv'],surrogates_subsample)
        end
    else
        if exist([surf_temp '/surrogates_subsample_PD_' hemi '.csv'], 'file')
            surrogates_subsample=csvread([surf_temp '/surrogates_subsample_PD_' hemi '.csv']);
        else
            surrogates_subsample_hipp = obj_subsample_hipp.fit(zscore(surf.coord(2,1:419))', n_surro);
            surrogates_subsample_DG = obj_subsample_DG.fit(zscore(surf.coord(2,420:483))', n_surro);
            surrogates_subsample=[surrogates_subsample_hipp;surrogates_subsample_DG];
            csvwrite([surf_temp '/surrogates_subsample_PD_' hemi '.csv'],surrogates_subsample)
        end    
    end   
    
    if g==3
        tbl=table(gradient_avg(:,g),zscore(surf.coord(1,:))','VariableNames',{'G3','AP'});
        lm=fitlm(tbl,'G3~1+AP+AP^2');
        Rsq_real = lm.Rsquared.Adjusted;
        r2_all(g,h)=Rsq_real;
        Rsq_surrogate=nan(n_surro,1);
        parfor i=1:n_surro            
            tbl=table(gradient_avg(:,g),surrogates_subsample(:,i),'VariableNames',{'G3','AP'});
            lm=fitlm(tbl,'G3~1+AP+AP^2');
            Rsq_surrogate(i) = lm.Rsquared.Adjusted;        
        end
        prctile_rank = mean(Rsq_real > Rsq_surrogate);
    else
        r_real = corr(zscore(surf.coord(2,:))',gradient_avg(:,g));
        r_surrogate = corr(surrogates_subsample, gradient_avg(:,g));
        prctile_rank = mean(r_real > r_surrogate);
        r2_all(g,h)=r_real.^2;  
    end    
    
   if prctile_rank>0.5
      p_all(g,h)=1-prctile_rank;
   else
      p_all(g,h)=prctile_rank;
   end
   
%    if g==1 || g==2 
%        f=figure;h1=histfit(r_surrogate,100);set(h1(1),'facecolor','w');
%        hold on;plot([r_real,r_real],[0,20],'y-','LineWidth',4)
%        title(['G' num2str(g) '(' hemi ')- PD coord: r=' num2str(r_real) ', p=' num2str(p_all(g,h))])
%        saveas(f, [G_scatter_figure '/G' num2str(g) '(' hemi ')_PD_coord'],'tif')    
%    else 
%        f=figure;h1=histfit(Rsq_surrogate,100);set(h1(1),'facecolor','w');
%        hold on;plot([Rsq_real,Rsq_real],[0,20],'y-','LineWidth',4)       
%        title(['G3(' hemi ')- AP coord: R-sq-adj=' num2str(Rsq_real) ', p=' num2str(p_all(g,h))])
%        saveas(f, [G_scatter_figure '/G' num2str(g) '(' hemi ')_AP_coord'],'tif')      
%    end
   
end
end

clear f
for h=1:2
    hemi=hemis{h};
    surf=SurfStatReadSurf([surf_temp '/tpl-avg_' hemi '_space-unfold_den-2mm_midthickness.obj']);   
    gradient_avg=-gradient_2{h};
for g=1:3 

    if g==1
        f(1,2*g+h-2)=gramm('x',zscore(surf.coord(2,:)), 'y',gradient_avg(:,g));
        f(1,2*g+h-2).set_names('x','z-scored P-D coordinate','y',['G' num2str(g) '-' hemi '-grp'])
        f(1,2*g+h-2).geom_point(); f(1,2*g+h-2).stat_glm(); 
        f(1,2*g+h-2).set_point_options('markers',{'o'});
        if p_all(g,h)==0 
          f(1,2*g+h-2).set_title(['Rsq=' num2str(r2_all(g,h)) ', p<0.001']) 
        else
          f(1,2*g+h-2).set_title(['Rsq=' num2str(r2_all(g,h)) ', p=' num2str(p_all(g,h))])
        end
    elseif g==2
        f(1,2*g+h-2)=gramm('x',zscore(surf.coord(2,:)), 'y',gradient_avg(:,g));
        f(1,2*g+h-2).set_names('x','z-scored P-D coordinate','y',['G' num2str(g) '-' hemi '-grp'])
        f(1,2*g+h-2).geom_point(); f(1,2*g+h-2).stat_glm();
        f(1,2*g+h-2).set_point_options('markers',{'o'});
        if p_all(g,h)==0 
          f(1,2*g+h-2).set_title(['Rsq=' num2str(r2_all(g,h)) ', p<0.001']) 
        else
          f(1,2*g+h-2).set_title(['Rsq=' num2str(r2_all(g,h)) ', p=' num2str(p_all(g,h))])
        end
    else
        f(1,2*g+h-2)=gramm('x',zscore(surf.coord(1,:)), 'y',gradient_avg(:,g));
        f(1,2*g+h-2).set_names('x','z-scored A-P coordinate','y',['G' num2str(g) '-' hemi '-grp'])
        f(1,2*g+h-2).geom_point(); f(1,2*g+h-2).stat_fit('fun',@(b,c,d,x)b*x.^2+c*x+d,'intopt','functional');
        f(1,2*g+h-2).set_point_options('markers',{'o'});
        if p_all(g,h)==0 
          f(1,2*g+h-2).set_title(['Rsq=' num2str(r2_all(g,h)) ', p<0.001']) 
        else
          f(1,2*g+h-2).set_title(['Rsq=' num2str(r2_all(g,h)) ', p=' num2str(p_all(g,h))])
        end
    end      
    
end
end
f.set_title(['Visualization of the relationship between Gradient and AP(or PD) coordinate'])
f.set_color_options('map','brewer2')
fig=figure('Position',[100,100,1000,200]);
f.draw();
f.export('file_name','G-corr-AP_PD_suppl','export_path',[G_scatter_figure '/'],'file_type','pdf','width',60,'height',10)

%%
function [] = parsave(dir,conn_mat)
% save x in dir
% so I can save in parfor loop
save(dir,'conn_mat');
end