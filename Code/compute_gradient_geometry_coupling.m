clc;clear;
cd E:\code\hippocampus_cortex_gradient_youth % replace with absolute path of your work directory

addpath(genpath('./Dependencies/Matlab/'))

result_dir=['./Results/GeoEigenMode/'];
if ~exist(result_dir, 'dir')
    mkdir(result_dir)
end

% Path to the geometric eigenmodes
G_indiv_dir='./Data/hippovol_native/';

hemis={'L','R'};

opts=detectImportOptions('./code/HCD_LS_2.0_subject_completeness_bh.csv');
demographic_info=readtable('./code/HCD_LS_2.0_subject_completeness_bh.csv',opts);
subjlist=demographic_info{:,1};
subj_num=length(subjlist);

% initialize a cell to store the coupling coefficients between hippocampal
% gradients and geometric eigenmodes
feat_Eig_G={nan(subj_num,3),nan(subj_num,3)};

spaces={'canonical', 'unfold'};
surf_temp='./Surf_temp';

num_vert=483; num_mode=31;

surf=SurfStatReadSurf([surf_temp '/tpl-avg_L_space-unfold_den-2mm_midthickness.obj']);

if ~exist([result_dir '/feat_Eig_G.mat'],'file')
    for h=1:2
        hemi=hemis{h};
        gradient_avg=importdata(['./Results/gradient_avg/avg-gradient-4-' hemi '-652-dm.mat']);
        gradient_avg=-gradient_avg;
        eigenmode_sum=zeros(num_vert,num_mode);
        k=0;
        for i=1:subj_num
            subj_id=subjlist{i};

            %[G_indiv_dir '/sub-' subj_id '/sub-' subj_id '_hemi-' hemi '_space-T1w_desc-hippo-2mm_emode_31.txt']
            %[G_indiv_dir '/sub-' subj_id '/sub-' subj_id '_hemi-' hemi ...
        %         '_space-T1w_den-2mm_label-hipp_midthickness_emode_200.txt']
            hipp_mode_file= [G_indiv_dir '/sub-' subj_id '/sub-' subj_id '_hemi-' hemi ...
                '_space-T1w_desc-hipp-2mm_emode.func.gii'];
            DG_mode_file= [G_indiv_dir '/sub-' subj_id '/sub-' subj_id '_hemi-' hemi ...
                '_space-T1w_desc-dentate-2mm_emode.func.gii'];    

            if ~exist(hipp_mode_file,'file')
               fprintf([subj_id,'-',hemi,'-hipp not exists\n']);
               hipp_mode=gifti([[G_indiv_dir '/sub-' subj_id '/sub-' subj_id '_hemi-' hemis{3-h} ...
                '_space-T1w_desc-hipp-2mm_emode.func.gii']]);
               hipp_mode=hipp_mode.cdata;
        %        continue;
            else
               hipp_mode=gifti(hipp_mode_file);hipp_mode=hipp_mode.cdata; 
            end
        %     hipp_mode=gifti(hipp_mode_file);hipp_mode=hipp_mode.cdata;

            if ~exist(DG_mode_file,'file')
        %         dg_mode=nan(64,31);
                dg_mode=gifti([G_indiv_dir '/sub-' subj_id '/sub-' subj_id '_hemi-' hemis{3-h} ...
                '_space-T1w_desc-dentate-2mm_emode.func.gii']);
                dg_mode=dg_mode.cdata;
                fprintf([subj_id,'-',hemi,'-DG not exists\n']);
            else
                dg_mode=gifti(DG_mode_file);dg_mode=dg_mode.cdata;
            end

            eigenmode=[hipp_mode;dg_mode];
            for j=1:10
                ind=~(isnan(eigenmode(:,j+1)) | isnan(gradient_avg(:,j)));
                if corr(eigenmode(ind,j+1),gradient_avg(ind,j)) <0
                    eigenmode(ind,j+1)=-eigenmode(ind,j+1);
                end
            end

            if isempty(find(isnan(eigenmode), 1))
                eigenmode_sum=eigenmode_sum+eigenmode;
                k=k+1;
            end

            M1=eigenmode(:,2)'; save([result_dir '/Individual/M1/' subj_id '_' hemi '_M1.mat'],'M1');
            M2=eigenmode(:,3)'; save([result_dir '/Individual/M2/' subj_id '_' hemi '_M2.mat'],'M2');
            M3=eigenmode(:,4)'; save([result_dir '/Individual/M3/' subj_id '_' hemi '_M3.mat'],'M3');

            gradient_indiv=importdata(['./Results/gradient_indiv/gradient_indiv(sm0)/' subj_id '_V1_MR-' hemi '-gradient.mat']);
            gradient_indiv=-gradient_indiv;

            feat_Eig_G{h}(i,1)=corr(M1',gradient_indiv(:,1));feat_Eig_G{h}(i,2)=corr(M2',gradient_indiv(:,2));
            feat_Eig_G{h}(i,3)=corr(M3',gradient_indiv(:,3));

        end

        eigenmode_avg=eigenmode_sum/k;
        if ~exist([result_dir, '/eigenmode_avg_' hemi '.csv'],'file')    
            csvwrite([result_dir, '/eigenmode_avg_' hemi '.csv'],eigenmode_avg)
        end

        % if ~exist([result_dir '/Figs/Eigenmode_corr_G' '-', hemi '.tif'], 'file')
        %     G_corr_mode=corr(gradient_avg,eigenmode_avg);
        %     csvwrite([result_dir, '/eigenmode_corr_gradient_avg_' hemi '.csv'],G_corr_mode)
        % 
        %     f=figure;ax=imagesc(abs(G_corr_mode(:,2:11)),[0 1]);xticks(1:1:10)
        %     xlabel('Eigenmode');ylabel('gradient');colormap(bluewhitered);colorbar
        %     saveas(f, [result_dir '/Figs/Eigenmode_corr_G' '-', hemi],'tif')
        % end
    end
    save([result_dir '/feat_Eig_G.mat'], 'feat_Eig_G')
end


if ~exist([result_dir '/Figs'],'dir')
    mkdir([result_dir '/Figs'])
end
%% [Figure 1g] Visualize the averaged hippocampal geometric eigenmode
for h=1:2
    hemi=hemis{h};
    eigenmode_avg=csvread([result_dir, '/eigenmode_avg_' hemi '.csv']);
for sp=1:2
    space=spaces{sp};
    surf=SurfStatReadSurf([surf_temp '/tpl-avg_space-' space '_den-2mm_midthickness.obj']);
    if sp==1 && h==1
       surf.coord(1,:)=max(surf.coord(1,:))+0.5-surf.coord(1,:);
    elseif sp==2 && h==2
       surf.coord(2,:)=max(surf.coord(2,:))+0.5-surf.coord(2,:);
    end
    
    for m=2:4
        figure('Color','White','Units','Normalized','Position',[.1 .1 .56 .66]); 
        [a,~]=BladeSurfStatViewData(eigenmode_avg(:,m)', surf, ['Eigenmode_', num2str(m-1), '-', hemi, '_', space, '-avg']);
        colormap(bluewhitered);camlight;SurfStatColLim([min(eigenmode_avg(:,m)) max(eigenmode_avg(:,m))]);
        colorbar('ytick',linspace(min(eigenmode_avg(:,m)), max(eigenmode_avg(:,m)),6));
        if sp==2;view(-90,90);end
        saveas(a, [result_dir '/Figs/Eigenmode_', num2str(m-1), '-', hemi, '_', space],'tif')
%         colorbar off
%         saveas(a, [result_dir '/Figs/Eigenmode_', num2str(m-1), '-', hemi, '_', space],'pdf')
    end    
end
end

%% [Figure 1h] Test the relationships between hipocampal gradients and geometric eigenmodes
geo_dist_hipp=gifti([surf_temp '/hipp/midthickness_geoDist.gii']);geo_dist_hipp=geo_dist_hipp.cdata;
geo_dist_DG=gifti([surf_temp '/dentate/midthickness_geoDist.gii']);geo_dist_DG=geo_dist_DG.cdata;

rand_init=0;n_surro=1000;
num_samp1=400;
num_neigh1=80;
num_samp2=60;
num_neigh2=10;

obj_subsample_hipp = variogram(geo_dist_hipp, 'ns', num_samp1, ...
    'knn', num_neigh1, 'random_state', rand_init, 'num_workers', 48);
obj_subsample_DG = variogram(geo_dist_DG, 'ns', num_samp2, ...
    'knn', num_neigh2, 'random_state', rand_init, 'num_workers', 48);
 
p_all=nan(3,2);
r_all=nan(3,2);

for h=1:2
 hemi=hemis{h};
 gradient_avg=importdata(['./Results/gradient_avg/avg-gradient-4-' hemi '-652-dm.mat']);
 gradient_avg=-gradient_avg; 
 
 eigenmode_avg=importdata([result_dir, '/eigenmode_avg_' hemi '.csv']);
 
    for g=1:3
        if exist([surf_temp '/surrogates_subsample_mode' num2str(g+1) '_' hemi '.csv'], 'file')
            surrogates_subsample=csvread([surf_temp '/surrogates_subsample_mode' num2str(g+1) '_' hemi '.csv']);
        else        
            surrogates_subsample_hipp = obj_subsample_hipp.fit(eigenmode_avg(1:419,g+1), n_surro);
            surrogates_subsample_DG = obj_subsample_DG.fit(eigenmode_avg(420:483,g+1), n_surro);
            surrogates_subsample=[surrogates_subsample_hipp;surrogates_subsample_DG];
            csvwrite([surf_temp '/surrogates_subsample_mode' num2str(g+1) '_' hemi '.csv'],surrogates_subsample)
        end
        r_real = corr(eigenmode_avg(:,g+1),gradient_avg(:,g));
        r_surrogate = corr(surrogates_subsample, gradient_avg(:,g));
        prctile_rank = mean(r_real > r_surrogate);
        r_all(g,h)=r_real;
        
        if prctile_rank>0.5
           p_all(g,h)=1-prctile_rank;
        else
           p_all(g,h)=prctile_rank;
        end
        
%        fig=figure;h1=histfit(r_surrogate,100);set(h1(1),'facecolor','w');
%        hold on;plot([r_real,r_real],[0,20],'y-','LineWidth',4)
%        title(['G' num2str(g) '(' hemi ') - EigenMode: r=' num2str(r_real) ', p=' num2str(p_all(g,h))])
%        saveas(fig, [result_dir '/Figs/G' num2str(g) '(' hemi ')_EigenMode'],'tif') 
    end    
end

clear f
for h=1:2
 hemi=hemis{h};
 gradient_avg=importdata(['./Results/gradient_avg/avg-gradient-4-' hemi '-652-dm.mat']);
 gradient_avg=-gradient_avg; 
 
 eigenmode_avg=importdata([result_dir, '/eigenmode_avg_' hemi '.csv']);
for g=1:3 

    f(1,2*g+h-2)=gramm('x',eigenmode_avg(:,g+1), 'y',gradient_avg(:,g));
    f(1,2*g+h-2).set_names('x',['Geometric eigenmode' num2str(g)],'y',['G' num2str(g) '-' hemi '-grp'])
    f(1,2*g+h-2).geom_point();f(1,2*g+h-2).stat_glm(); 
    if p_all(g,h)==0 
      f(1,2*g+h-2).set_title(['r=' num2str(r_all(g,h)) ', p<0.001']) 
    else
      f(1,2*g+h-2).set_title(['r=' num2str(r_all(g,h)) ', p=' num2str(p_all(g,h))])
    end
end
end
f.set_title(['Visualization of the relationship between Gradient and Geometric eigenmode'])
f.set_color_options('map','brewer2')
fig=figure('Position',[100,100,1200,800]);
f.draw();
f.export('file_name','G-corr-mode','export_path',[result_dir '/Figs/'],'file_type','pdf','width',60,'height',10)