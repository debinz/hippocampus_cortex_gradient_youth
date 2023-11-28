clc;clear;
cd E:\code\hippocampus_cortex_gradient_youth %replace with absolute path of your work directory

addpath(genpath('./Dependencies/Matlab/'))

is_conn_thresh=false; %true;%
conn_sparsity=90;
hemis={'L','R'};
vert_num=360;  % number of cortical regions/vertices

avg_result='./Results/gradient_avg';
g_indiv='./Results/gradient_indiv/';
conn_indiv='./Results/conn_mat/';

surf_temp = './Surf_temp';

root_dir='./Results/hipp_CTX_proj';
if ~exist(root_dir, 'dir')
    mkdir(root_dir)
end

result_folder=[root_dir '/indiv' num2str(vert_num)];

if ~exist(result_folder,'dir')
   mkdir(result_folder) ;
end

fig_dir=[root_dir '/figure' num2str(vert_num)];    
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

%% compute the gradient projections
% read the demographic information of HCP-D data
opts=detectImportOptions('./Code/HCD_LS_2.0_subject_completeness_bh.csv');
demographic_info=readtable('./Code/HCD_LS_2.0_subject_completeness_bh.csv',opts);
subjlist=demographic_info{:,1};
subj_num=length(subjlist);

% load cortical functional hierarchy map
g1_360=importdata([surf_temp '/cortex_surf/g1_360.mat']);

G_CTX_grp_LR=cell(1,2);
G1_CTX_all_LR=cell(1,2);G2_CTX_all_LR=cell(1,2);G3_CTX_all_LR=cell(1,2);
for h=1:2
    hemi=hemis{h};
    
    if ~exist([result_folder,'/G_CTX_grp_' hemi '.csv'],'file')
    
    k=0;
    G1_CTX_all=nan(subj_num,vert_num);G2_CTX_all=nan(subj_num,vert_num);G3_CTX_all=nan(subj_num,vert_num);
    G_CTX=zeros(vert_num,10);
    G1_all=nan(subj_num,483);G2_all=nan(subj_num,483);G3_all=nan(subj_num,483);
    G_CTX_FH=nan(subj_num,6);
    parfor i=1:subj_num
        subj_id=[subjlist{i} '_V1_MR'];
        if ~exist([g_indiv '/' subj_id '-' hemi '-gradient.mat'],'file')
            disp([subj_id,'- lacking gradient'])
            continue
        else
            disp(['processing ' subj_id])
            if ~exist([result_folder '/' subj_id  '_' hemi '_CTX.mat'], 'file')
                C_r=importdata([conn_indiv '/' subj_id '_' hemi '_conn_mat.mat']);
                C_r=C_r';
                if is_conn_thresh
                   C_r(C_r < prctile(C_r,conn_sparsity)) = 0;  
                end
                gradient=importdata([g_indiv '/' subj_id '-' hemi '-gradient.mat']);            
                gradient=-gradient; % inversed the sign
                G1_all(i,:)=gradient(:,1);G2_all(i,:)=gradient(:,2);G3_all(i,:)=gradient(:,3);
                G_CTX_sbj=C_r*gradient;
                parsave([result_folder '/' subj_id  '_' hemi '_CTX.mat'],G_CTX_sbj);                                
            else
                G_CTX_sbj=importdata([result_folder '/' subj_id  '_' hemi '_CTX.mat']);
            end
            
            G_CTX_FH(i,:)=[corr(G_CTX_sbj(:,1),g1_360),...
                    corr(G_CTX_sbj(:,2),g1_360),...
                    corr(G_CTX_sbj(:,3),g1_360),...
                    corr(G_CTX_sbj(:,1),g1_360,'type','Spearman'),...
                    corr(G_CTX_sbj(:,2),g1_360,'type','Spearman'),...
                    corr(G_CTX_sbj(:,3),g1_360,'type','Spearman')];      
            
            G1_CTX_all(i,:)=G_CTX_sbj(:,1);G2_CTX_all(i,:)=G_CTX_sbj(:,2);G3_CTX_all(i,:)=G_CTX_sbj(:,3);
            G_CTX=G_CTX+G_CTX_sbj;
            k=k+1;            
        end                            

    end
    G_CTX_grp=G_CTX/k;
    G_CTX_grp_LR{h}=G_CTX_grp;
    csvwrite([result_folder,'/G_CTX_grp_' hemi '.csv'],G_CTX_grp)
    csvwrite([result_folder,'/G1_CTX_all_' hemi '.csv'],G1_CTX_all)
    csvwrite([result_folder,'/G2_CTX_all_' hemi '.csv'],G2_CTX_all)
    csvwrite([result_folder,'/G3_CTX_all_' hemi '.csv'],G3_CTX_all)
    csvwrite([g_indiv, '/G1_all_' hemi '.csv'],G1_all)
    csvwrite([g_indiv, '/G2_all_' hemi '.csv'],G2_all)
    csvwrite([g_indiv, '/G3_all_' hemi '.csv'],G3_all)
    csvwrite([result_folder, '/G_CTX_FH_all_' hemi '.csv'],G_CTX_FH)
    
    else
     G_CTX_grp=csvread([result_folder,'/G_CTX_grp_' hemi '.csv']);   
     G_CTX_grp_LR{h}=G_CTX_grp;
     G1_CTX_all=csvread([result_folder,'/G1_CTX_all_' hemi '.csv']);
     G1_CTX_all_LR{h}=G1_CTX_all;
     G2_CTX_all=csvread([result_folder,'/G2_CTX_all_' hemi '.csv']);
     G2_CTX_all_LR{h}=G2_CTX_all;
     G3_CTX_all=csvread([result_folder,'/G3_CTX_all_' hemi '.csv']);
     G3_CTX_all_LR{h}=G3_CTX_all;     
    end    
end
G_CTX_all_LR={G1_CTX_all_LR,G2_CTX_all_LR,G3_CTX_all_LR};

%% [Figure 1d] plot the averaged cortical distribution of hippocampal gradient projections
surf_ctx=SurfStatReadSurf([surf_temp '/cortex_surf/surf_cortex_temp_inflated_LR.obj']);
[surf_lh,surf_rh]=split_surfaces(surf_ctx);
% [surf_lh,surf_rh]=load_conte69();
labeling = load_parcellation('glasser',360);
for h=1:2
    hemi=hemis{h};
    
    G_CTX_grp=G_CTX_grp_LR{h};
       
    f=plot_hemispheres(G_CTX_grp(:,1),{surf_lh,surf_rh},'parcellation',labeling.glasser_360);
    f.colormaps([0.8,0.8,0.8;viridis])
    f.colorlimits([min(G_CTX_grp(:,1))*1.01, max(G_CTX_grp(:,1))]);
    saveas(f.handles.figure, [fig_dir '/G1(' hemi ')_CTX_grp'],'tif')
    
    f=plot_hemispheres(G_CTX_grp(:,2),{surf_lh,surf_rh},'parcellation',labeling.glasser_360);
    f.colormaps([0.8,0.8,0.8;viridis])
    f.colorlimits([min(G_CTX_grp(:,2))*1.02, max(G_CTX_grp(:,2))]);
    saveas(f.handles.figure, [fig_dir '/G2(' hemi ')_CTX_grp'],'tif')
    
    f=plot_hemispheres(G_CTX_grp(:,3),{surf_lh,surf_rh},'parcellation',labeling.glasser_360);
    f.colormaps([0.8,0.8,0.8;viridis])
    f.colorlimits([min(G_CTX_grp(:,3))*1.001, max(G_CTX_grp(:,3))]);
    saveas(f.handles.figure, [fig_dir '/G3(' hemi ')_CTX_grp'],'tif')
end

%% [Figure 1e] Select top and bottom 5% projection value to show the pattern of cortical connectivity along the hippocampal gradient
% make the color map for this figure
cm=viridis(360);
ind=[round(size(cm,1)*0.05),round(size(cm,1)*0.95)];
cm(ind(1)+1:ind(2)-1,:)=repmat([0.8,0.8,0.8],ind(2)-ind(1)-1,1);
% figure;rgbplot(cm(1:2*ind(1),:));hold on;colormap(cm(1:2*ind(1),:));colorbar('Ticks',[])
% figure;rgbplot(cm(ind(2)-ind(1):end,:));hold on;colormap(cm(ind(2)-ind(1):end,:));colorbar('Ticks',[])

clear f
figure('Position',[100,100,1500,800]);

for h=1:2    
    hemi=hemis{h};
    conn_mean=importdata([avg_result, '/conn_matrix_mean_' hemi '-652.mat']);
    g_avg=importdata([avg_result, '/avg-gradient-4-' hemi '-652-dm.mat']);
    G_CTX_grp=csvread([result_folder,'/G_CTX_grp_' hemi '.csv']);

    if is_conn_thresh
        conn_mean(conn_mean < prctile(conn_mean,conn_sparsity,2)) = 0;  
    end

    for g=1:3
        [g_sort,ind]=sort(-g_avg(:,g)); g_sort=g_sort';
        conn_sortbyy=conn_mean(ind,:); conn_y_all=mat2cell(conn_sortbyy',ones(360,1),483);
        prj_g=G_CTX_grp(:,g); prj_g_cell=num2cell(prj_g)';
    %     roi_ind=[1:1:360];     

        f(h,g)=gramm('x',g_sort,'y',conn_y_all,'color',prj_g, ...
            'subset', prj_g< prctile(prj_g,5) ...
            |  prj_g>prctile(prj_g,95));
        f(h,g).geom_point();
    %     f(h,g).stat_smooth();

        f(h,g).set_names('x',['Hipp_Gradient ' num2str(g)], 'y', 'Cortical connectivity', 'color', ['G' num2str(g) '_proj'])
        f(h,g).set_color_options('map',cm) %viridis
    end
end
f.set_title(['Cortical connectivity along the hippocampal gradients (select top and bottom 5% projection value)'])
f.draw();
f.export('file_name','CTX_conn_grp_G','export_path',fig_dir,'file_type','pdf','width',45,'height',20)
f.export('file_name','CTX_conn_grp_G','export_path',fig_dir,'file_type','png','width',45,'height',20)

%% [Figure 1f] study the relationships between hippocampal gradient projections and cortical functional hierarchy
geo_dist_L=importdata([surf_temp '/cortex_surf/LeftParcelGeodesicDist.mat']);
geo_dist_R=importdata([surf_temp '/cortex_surf/RightParcelGeodesicDist.mat']);

rand_init=0;n_surro=1000;
num_samp=160;
num_neigh=80;

% for both two cortical hemis
p_all=nan(3,2);
mkdir([fig_dir '/permutation-FH/']);

if ~exist([surf_temp '/cortex_surf/surrogates_subsample_CTX-g1_360.csv'],'file')
    obj_subsample_L = variogram(geo_dist_L, 'ns', num_samp, ...
        'knn', num_neigh, 'random_state', rand_init);
    surrogates_subsample_L = obj_subsample_L.fit(g1_360(1:180,1), n_surro);

    obj_subsample_R = variogram(geo_dist_R, 'ns', num_samp, ...
        'knn', num_neigh, 'random_state', rand_init);
    surrogates_subsample_R = obj_subsample_R.fit(g1_360(181:360,1), n_surro);

    surrogates_subsample=[surrogates_subsample_L;surrogates_subsample_R];
    csvwrite([surf_temp '/cortex_surf/surrogates_subsample_CTX-g1_360.csv'],surrogates_subsample);
else
    surrogates_subsample=csvread([surf_temp '/cortex_surf/surrogates_subsample_CTX-g1_360.csv']);
end

for h=1:2
    hemi=hemis{h}; 
    
for g=1:3    

    r_real = corr(g1_360(:,1), G_CTX_grp_LR{h}(:,g));
    r_surrogate = corr(surrogates_subsample, G_CTX_grp_LR{h}(:,g));
    prctile_rank = mean(r_real > r_surrogate);
    
    if prctile_rank>0.5
      p_all(g,h)=1-prctile_rank;
    else
      p_all(g,h)=prctile_rank;
    end

%     f=figure;h1=histfit(r_surrogate,100);set(h1(1),'facecolor','w');
%     hold on;plot([r_real,r_real],[0,20],'y-','LineWidth',4)
%     title(['hipp(' hemi ')- CTX projection: r=' num2str(r_real) ', p=' num2str(prctile_rank)])
%     saveas(f, [fig_dir '/permutation-FH/G' num2str(g) '(' hemi ')_CTX_FH'],'tif')    

end
end

clear f 
for h=1:2
    hemi=hemis{h};
    for g=1:3
               
       f(g,h)=gramm('x',g1_360(:,1), 'y',G_CTX_grp_LR{h}(:,g));
       f(g,h).geom_point();f(g,h).stat_glm();
       f(g,h).set_names('x','FH','y',['G' num2str(g) '-' hemi '-CTX-projection'])
       [r,~]=corr(g1_360(:,1),G_CTX_grp_LR{h}(:,g));
       if p_all(g,h)==0 
          f(g,h).set_title(['r=' num2str(r) ', p<0.001']) 
       else
          f(g,h).set_title(['r=' num2str(r) ', p=' num2str(p_all(g,h))])
       end
    end
end
f.set_title(['Relationship between Gradient projection and Functional hierarchy'])
f.set_color_options('map','brewer2')
fig=figure('Position',[100,100,800,1200]);
f.draw();
f.export('file_name','G-corr-FH','export_path',[fig_dir '/permutation-FH/'],'file_type','pdf','width',20,'height',30)
f.export('file_name','G-corr-FH','export_path',[fig_dir '/permutation-FH/'],'file_type','png','width',20,'height',30)

%%
function [] = parsave(dir,G_CTX_sbj)
% save x in dir
% so I can save in parfor loop
save(dir,'G_CTX_sbj');
end
