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

age_grps=[5 14 22]; 

for grp=1:(length(age_grps)-1)
        
    %%get the behav residuals that controled the sex effect
    dccs_res=readtable([bh_dir '/Fig_dev/residual/dccs_residuals_age.csv']);dccs_res=table2array(dccs_res);
    flanker_res=readtable([bh_dir '/Fig_dev/residual/flanker_residuals_age.csv']);flanker_res=table2array(flanker_res);
    wm_res=readtable([bh_dir '/Fig_dev/residual/wm_residuals_age.csv']);wm_res=table2array(wm_res);

    bh_res={dccs_res,flanker_res,wm_res};
    bh_var={'dccs','flanker','wm'};

    %get the s-f coupling residuals that controled the sex, site, and head motion(mFD) effect
    EigG1_res=readtable([sf_dir '/Fig_dev/residual/EigMod_G1_residuals_feat_age.csv']);EigG1_res=table2array(EigG1_res);
    EigG2_res=readtable([sf_dir '/Fig_dev/residual/EigMod_G2_residuals_feat_age.csv']);EigG2_res=table2array(EigG2_res);
    EigG3_res=readtable([sf_dir '/Fig_dev/residual/EigMod_G3_residuals_feat_age.csv']);EigG3_res=table2array(EigG3_res);
    EigG_res={EigG1_res,EigG2_res,EigG3_res};

    r_all=nan(3,6);p_all=nan(3,6);sig_all=nan(3,6);
    hemis={'L','R'};
    for b=1:3
        ind=demographic_info.Age>age_grps(grp) & demographic_info.Age<age_grps(grp+1) & ~isnan(demographic_info.(bh_var{b}));
        ind_bh=demographic_info.Age(~isnan(demographic_info.(bh_var{b})))>age_grps(grp) & ...
            demographic_info.Age(~isnan(demographic_info.(bh_var{b})))<age_grps(grp+1);
        for g=2:2
            for h=1:2
                [r_all(b,2*g-2+h),p_all(b,2*g-2+h),sig_all(b,2*g-2+h),fig]= ...
                  permutation_test(bh_res{b}(ind_bh,1),EigG_res{g}(ind,h),...
                     1000,{bh_var{b},['EigG' num2str(g) '-' hemis{h}]},'Pearson'); % 'Spearman''Pearson''Kendall' 
%                 saveas(fig,[bh_dir '/Fig_BrBeAssoc/age_grp/' bh_var{b} '-' 'EigG' num2str(g) '_' hemis{h} '_' num2str(age_grps(grp))],'tif')  
            end
        end    
    end

    clear f
    fig=figure('position',[200,200,800,300]);
    for b=1:3
        ind=demographic_info.Age>age_grps(grp) & demographic_info.Age<age_grps(grp+1) & ~isnan(demographic_info.(bh_var{b}));
        ind_bh=demographic_info.Age(~isnan(demographic_info.(bh_var{b})))>age_grps(grp) & ...
            demographic_info.Age(~isnan(demographic_info.(bh_var{b})))<age_grps(grp+1); 
        for g=2:2
            for h=1:2
                f(h,b)=gramm('x',EigG_res{g}(ind,h),'y',bh_res{b}(ind_bh,1));
                f(h,b).geom_point();f(h,b).stat_glm();
                f(h,b).set_names('x',['EigG' num2str(g) '-' hemis{h}],'y',bh_var{b});
                if p_all(b,2*g-2+h)==0
                    f(h,b).set_title(['r = ' num2str(r_all(b,2*g-2+h)) ', p < 0.001 '])
                else
                    f(h,b).set_title(['r = ' num2str(r_all(b,2*g-2+h)) ', p = ' num2str(p_all(b,2*g-2+h))])
                end
            end
        end
    end
    f.set_title(['Relationships between Eigenmod-Gradient2 couping and Excutive Function performance_age' num2str(age_grps(grp)) '~' num2str(age_grps(grp+1)-1)])
    f.set_color_options('map','brewer2')
    f.draw();
    f.export('file_name',['sfcouping(g' num2str(g) ')_bh_' num2str(age_grps(grp))],'export_path',Fig_bh,'file_type','pdf','width',30,'height',20);
%     saveas(fig,[bh_dir '/Fig_BrBeAssoc/age_grp/sfcouping(g' num2str(g) ')_bh_' num2str(age_grps(grp))],'tif')  


end

function [r,p_value,sig,fig]=permutation_test(x,y,nIterNull,var_names,type)
    r=corr(x, y, 'type',type);
    r_rand_all=zeros(nIterNull,1);
    for i=1:nIterNull
        rand_ind=randperm(length(x));
        r_rand_all(i)=corr(x(rand_ind,1), y, 'type',type);
    end
    p_value=mean(r>r_rand_all);
    if p_value>0.5
        p_value=1-p_value;
    end
    sig=p_value<=0.025;
    fig=0;
%     fig=figure;h=histfit(r_rand_all,200);hold on;plot([r,r],[0,20],'y-','LineWidth',2);
%     xlabel([type ' correlations']);ylabel('count');set(h(1),'facecolor','w','EdgeColor','k');
%     title([var_names{1} '-' var_names{2} ': r=' num2str(r) ', p=' num2str(p_value)]);
%     box off
end

