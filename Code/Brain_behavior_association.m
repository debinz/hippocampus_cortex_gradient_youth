clc;clear;

Root_dir='E:\Github\hippocampus_cortex_gradient_youth'; % replace this with absolute path of your working directory

addpath(genpath([Root_dir '/Dependencies/Matlab/']))

bh_dir=[Root_dir '/Results/Cog_Behavior'];
Fig_bh=[bh_dir '/Fig_BrBhAssoc'];

if ~exist(Fig_bh, 'dir')
   mkdir(Fig_bh) 
end

sf_dir=[Root_dir '/Results/GeoEigenMode'];

% read the demographic information of HCP-D data
opts=detectImportOptions([Root_dir '/Code/HCD_LS_2.0_subject_completeness_bh.csv']);
demographic_info=readtable([Root_dir '/Code/HCD_LS_2.0_subject_completeness_bh.csv'],opts);
subjlist=demographic_info{:,1};
subj_num=length(subjlist);

%get the behav residuals that controled the sex effect
dccs_res=readtable([bh_dir '/Fig_dev/residual/dccs_residuals_age.csv']);dccs_res=table2array(dccs_res);
flanker_res=readtable([bh_dir '/Fig_dev/residual/flanker_residuals_age.csv']);flanker_res=table2array(flanker_res);
wm_res=readtable([bh_dir '/Fig_dev/residual/wm_residuals_age.csv']);wm_res=table2array(wm_res);

bh_res={dccs_res,flanker_res,wm_res};
bh_var={'dccs','flanker','wm'};

%get the s-f coupling residuals that controled the sex, scanning site, and head motion(mFD) effect
EigG1_res=readtable([sf_dir '/Fig_dev/residual/EigMod_G1_residuals_feat_age.csv']);EigG1_res=table2array(EigG1_res);
EigG2_res=readtable([sf_dir '/Fig_dev/residual/EigMod_G2_residuals_feat_age.csv']);EigG2_res=table2array(EigG2_res);
EigG3_res=readtable([sf_dir '/Fig_dev/residual/EigMod_G3_residuals_feat_age.csv']);EigG3_res=table2array(EigG3_res);
EigG_res={EigG1_res,EigG2_res,EigG3_res};

r_all=nan(3,6);p_all=nan(3,6);sig_all=nan(3,6);
hemis={'L','R'};
for b=1:3
    ind=~isnan(demographic_info.(bh_var{b}));
    for g=1:1
        for h=1:2
            [r_all(b,2*g-2+h),p_all(b,2*g-2+h),sig_all(b,2*g-2+h),fig]= ...
              permutation_test(bh_res{b}(:,1),EigG_res{g}(ind,h),...
                 1000,{bh_var{b},['EigG' num2str(g) '-' hemis{h}]});
%             saveas(fig,[bh_dir '/Fig_BrBeAssoc/' bh_var{b} '-' 'EigG' num2str(g) '_' hemis{h}],'tif')  
        end
    end    
end

clear f
fig=figure('position',[200,200,1200,200]);
for b=1:3
    ind=~isnan(demographic_info.(bh_var{b}));
    for g=1:1
        for h=1:2
            f(1,2*b+h-2)=gramm('x',EigG_res{g}(ind,h),'y',bh_res{b}(:,1));
            f(1,2*b+h-2).geom_point();f(1,2*b+h-2).stat_glm();
            f(1,2*b+h-2).set_names('x',['EigG1-' hemis{h}],'y',bh_var{b});
            if p_all(b,2*g-2+h)==0
                f(1,2*b+h-2).set_title(['r = ' num2str(r_all(b,2*g-2+h)) ', p < 0.001 '])
            else
                f(1,2*b+h-2).set_title(['r = ' num2str(r_all(b,2*g-2+h)) ', p = ' num2str(p_all(b,2*g-2+h))])
            end
        end
    end
end
f.set_title('Relationships between Eigenmod-Gradient couping and Excutive Function performance')
f.set_color_options('map','brewer2')
f.draw();
f.export('file_name','sf2couping_bh','export_path',Fig_bh,'file_type','pdf','width',60,'height',10);

function [r,p_value,sig,fig]=permutation_test(x,y,nIterNull,var_names)
    r=corr(x,y);
    r_rand_all=zeros(nIterNull,1);
    for i=1:nIterNull
        rand_ind=randperm(length(x));
        r_rand_all(i)=corr(x(rand_ind,1),y);
    end
    p_value=mean(r>r_rand_all);
    if p_value>0.5
        p_value=1-p_value;
    end
    sig=p_value<=0.025;
    fig=0;
%     fig=figure;h=histfit(r_rand_all,200);hold on;plot([r,r],[0,20],'y-','LineWidth',2);
%     xlabel('Pearson correlations');ylabel('count');set(h(1),'facecolor','w','EdgeColor','k');
%     title([var_names{1} '-' var_names{2} ': r=' num2str(r) ', p=' num2str(p_value)]);
%     box off
end
