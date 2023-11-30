clc;clear;
cd E:\hippo_gradient\HCD-unfold\
addpath(genpath('E:\code\'));
addpath('../')
% load('./figure360/cortex_surf/glasser360_sphere.mat');
% glasser360=csvread('./figure360/cortex_surf/glasser_360_conte69.csv');
surf_temp_LR=SurfStatReadSurf('./surf_temp/cortex_surf/surf_cortex_temp_inflated_LR.obj');
[surf_lh,surf_rh]=split_surfaces(surf_temp_LR);

result_folder='./G_CTX/figure360/ctx_dev_LR';
if ~exist(result_folder,'dir')
   mkdir(result_folder) 
end

% [surf_lh,surf_rh]=load_conte69();
labeling=load_parcellation('glasser',360);

yeo7_label=readtable('../glasser360_7networks.xlsx','ReadVariableNames',false);
% 
% f=plot_hemispheres(yeo7_label{:,3},{surf_lh,surf_rh},'parcellation',labeling.glasser_360);
% f.colormaps([0.8 0.8 0.8; viridis])
% f1.colorlimits([0.9 7]);
% g=1;
% hemi='L';
hemis={'L','R'};
yeo7=nan(7,3);

for g=1:3
%     for h=1:2
%         hemi=hemis{h};
%         
gam_results=readtable([result_folder '/G' num2str(g) '_CTX_proj_dev.csv']);
% age_deltaR2_sig=(str2double(gam_results.Anova_age_adjpvalue)<0.05).*gam_results.age_deltaR2;
range=minmax(gam_results.age_deltaR2(gam_results.age_adjpvalue<0.05)');
age_deltaR2_sig=(gam_results.age_adjpvalue<0.05).*gam_results.age_deltaR2;

for n=1:7
    ind=yeo7_label{yeo7_label{:,3}==n,1};
    yeo7(n,g)=mean(age_deltaR2_sig(ind));    
end

% if g==1
%     mm=[min(yeo7(:,g))*1.2, -min(yeo7(:,g))*1.2];
% elseif g==2
%     mm=[min(yeo7(:,g))*1.6, 0];
% else
%     mm=[min(yeo7(:,g))*1.2, -min(yeo7(:,g))*1.2];
% end

% if g==1
    mm=[0, max(yeo7(:,g))*1.2];
% elseif g==2
%     mm=[min(yeo7(:,g))*1.6, 0];
% else
%     mm=[min(yeo7(:,g))*1.2, -min(yeo7(:,g))*1.2];
% end

f=figure;

if g==1 || g==2
    De_mica_spider(yeo7(:,g),['Significant-deltaR2-Yeo7'],[mm(1), mm(2)],...
       {'Visual','Somatomotor','Dorsal Attention','Ventral Attention','Limbic','Frontoparietal','Default'},['G' num2str(g)],viridis,gca);
else
    De_mica_spider(yeo7(:,g)*10,['Significant-deltaR2-Yeo7'],[mm(1), mm(2)],...
       {'Visual','Somatomotor','Dorsal Attention','Ventral Attention','Limbic','Frontoparietal','Default'},['G' num2str(g)],viridis,gca);
end
% saveas(f, [result_folder '/Significant-deltaR2-Yeo7(G' num2str(g) ')'],'pdf')
   
age_deltaR2_sig(age_deltaR2_sig==0)=-2;
f1=plot_hemispheres(age_deltaR2_sig,{surf_lh,surf_rh},'parcellation',labeling.glasser_360);
f1.colormaps([0.8 0.8 0.8;viridis]) %flip(viridis)
% f1.colorlimits([range(1),range(2)]);
f1.colorlimits([0,max(gam_results.age_deltaR2)]);
% saveas(f1.handles.figure, [result_folder '/Delta_Rsq_sig_adj-hippo(G' num2str(g) ')'],'tif')

f2=plot_hemispheres(gam_results.age_deltaR2,{surf_lh,surf_rh},'parcellation',labeling.glasser_360);
f2.colormaps([0.8 0.8 0.8;viridis])
f2.colorlimits([min(gam_results.age_deltaR2)*200,max(gam_results.age_deltaR2)]);

saveas(f2.handles.figure, [result_folder '/Delta_Rsq_adj-hippo(G' num2str(g) ')'],'tif')
% 
% figure;plot(str2double(gam_results.Anova_age_pvalue), gam_results.age_pvalue,'.');
% xlabel('Anova.age.pvalue');ylabel('GAM.age.pvalue')
% hold on 
% % plot([0,0.1,1],[0,0.1,1],'y-')
% title(['G' num2str(g) ': relationship between gam.age.p and anova.age.p'])

% conn_64k_qvalue=zeros(64984,1);
% conn_64k_Rsq=zeros(64984,1);
% 
% for i=1:64984
%     if glasser360(i) ~= 0
%        conn_64k_qvalue(i)=gam_results.Anova_age_adjpvalue(glasser360(i)); 
%     else
%        conn_64k_qvalue(i)=20;
%     end
% end
% 
% for i=1:64984
%     if glasser360(i) ~= 0
%         if gam_results.Anova_age_adjpvalue(glasser360(i))<0.05
%            conn_64k_Rsq(i)=gam_results.GAM_age_deltaR2(glasser360(i)); 
%         else
%            conn_64k_Rsq(i)=20; 
%         end
%     else
%        conn_64k_Rsq(i)=20;
%     end
% end
% 
% deltaR2_range=minmax(gam_results.GAM_age_deltaR2');
% 
% f=figure;
% [a,~]=DeSurfStatViewData(conn_64k_Rsq,surf_temp_LR,['Delta-Rsq-adj: hippo(G' num2str(g) ')']);
% colormap([viridis;0.8 0.8 0.8]);SurfStatColLim([deltaR2_range(1),deltaR2_range(2)]);
% saveas(f, [result_folder '/Delta_Rsq_adj-hippo(G' num2str(g) ')'],'tif')

%     end
end

% figure;
% De_mica_spider(yeo7,['Significant-deltaR2-Yeo7'],[min(yeo7(:))*1.2, -min(yeo7(:))*1.2],...
%        {'Visual','Somatomotor','Dorsal Attention','Ventral Attention','Limbic','Frontoparietal','Default'},...
%        {'G1','G2','G3'},viridis(3),gca);