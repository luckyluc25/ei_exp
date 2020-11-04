%% Regress GABA/Gluatamte concentration on DV

%clean up
clear; close all; clc


%add the paths you need

%load dependent variables
cd( 'C:\Users\lucak\Desktop\clash_share-main\Second_Phase\DV');

files=dir('*.mat')
for k=1:numel(files)
    load(files(k).name);
end


%% read in GABA /glutamate concentrations
vp = {'1'; '2'; '3'; '4'; '5'; '7'; '9'; '10'; '12';  '14'; '15'; '16'; '17'; '18';  '20';  '22'; '23'; '24';  '26'; '27'; '28'; '29'; '30'; '31'; '32'; '33'; '34'; '35'; '38'};

% E/I balances are read out from csv files to a table and stored under % their MR name
table=importdata( 'C:\Users\lucak\Desktop\clash_share-main\ressources_share\table_anonym.mat');

%rename combines MR names with subject numbers
rename=importdata( 'C:\Users\lucak\Desktop\clash_share-main\ressources_share\rename_anonym.mat');

%in gm_spm12_test grey matter concentrations are stored by summing across %P(GM). They are later divided by total number of voxels in mask to
%approximate relative GM concentrations.
gm_spm12_gmwm=importdata('C:\Users\lucak\Desktop\clash_share-main\ressources_share\gm_spm12_gmwm_anonym.mat') ;

%just a check-if one variable is not there that should be, the code will
%break!
gm_spm12_gmwm.gm_all(gm_spm12_gmwm.gm_all==0)=NaN;
gm_spm12_gmwm.wm_all(gm_spm12_gmwm.wm_all==0)=NaN;
gm_spm12_gmwm.allvox(gm_spm12_gmwm.wm_all==0)=NaN;

%read in the following concentrations
labels_MRS={'LPFC','M1li','M1re','pgACC','aMCC'};


for concen=1:numel(labels_MRS)
    count=1;
    for it= 1:numel(vp)
        pbn=str2double(vp{it});
        if pbn<10
            subj_name=['clash0' num2str(pbn)];
        else
            subj_name=['clash' num2str(pbn)];
        end
        name=labels_MRS{concen};
        
        pbn_idx=(~cellfun(@isempty,(strfind(rename(:,5),subj_name))));
        mrName=rename(pbn_idx,6);mrName=regexp(mrName,'_','split');mrName=mrName{1,1};
        
        %find out where the mr name is in table, where the voxel name is
        %and whether this has been a valid measurement
        check1=(~cellfun(@isempty,(strfind(table.textdata.Alles(:,2),mrName))));
        check1(1)=[];
        check2=(~cellfun(@isempty,(strfind(table.textdata.Alles(:,2),name))));
        check2(1)=[];
        check3=table.data.Alles(:,1);
        
        idx=find(check1==1 & check2==1 & check3==1);
        
        
        %where is the correct GM concentration stored?
        spm12_idx=(~cellfun(@isempty,(strfind(gm_spm12_gmwm.mrName,mrName))));
        concen_con=find((~cellfun(@isempty,(strfind(erase(gm_spm12_gmwm.vox_name(spm12_idx,:),'_'),name)))));
        
        
        if isempty (idx)
            glut_corr(count,concen_con)=NaN;
            gaba_corr(count,concen_con)=NaN;
            count=count+1;
            continue
        end
        
        %Normalise by GM within each voxel
        
        glut_raw=table.data.Alles(idx,14);
        glut_raw=glut_raw/(gm_spm12_gmwm.gm_all(spm12_idx,concen_con)./gm_spm12_gmwm.allvox(spm12_idx,concen_con));
        gaba_raw=table.data.Alles(idx,10);
        gaba_raw=gaba_raw/(gm_spm12_gmwm.gm_all(spm12_idx,concen_con)./gm_spm12_gmwm.allvox(spm12_idx,concen_con));
        
        cr_raw=table.data.Alles(idx,8);
        cr_raw=cr_raw./((gm_spm12_gmwm.gm_all(spm12_idx,concen_con)+gm_spm12_gmwm.wm_all(spm12_idx,concen_con))/gm_spm12_gmwm.allvox(spm12_idx,concen_con));
        
        glut_corr(count,concen)=glut_raw/cr_raw;
        gaba_corr(count,concen)=gaba_raw/cr_raw;
        
        
        clear gaba_raw cr_raw glut_raw
        
        
        count=count+1;
        
        
    end
    
end

%estimate ratios
x=[(glut_corr./gaba_corr)];
x_single=[gaba_corr glut_corr];
model=importdata('C:\Users\lucak\Desktop\clash_final\Second_Phase\Modelling_lokal\extended_fitting\EVhybridlu.mat');
model_ddm=importdata( 'C:\Users\lucak\Desktop\clash_final\Second_Phase\DV\objective_allregs_vanda_short.mat');

depvar=model.params(:,5);

%============================================================================================================%
x = [ones(size(x,1),1) normalise([x sum_nb_obj ])];
x_single = [ones(size(x_single,1),1) normalise([x_single sum_nb_obj ])];
y = normalise(depvar);

lm = fitlm(x(:,2:size(x,2)),y)
coefCI(lm)

lm_single = fitlm(x_single(:,2:size(x_single,2)),y);


sig=find(lm.Coefficients.pValue<0.05);

%contrasts
tc=eye(size(x,2)); % matrix contrasts

% Are there any sigificant contributions- if not stop here.
if sum(sig)>0
    keyboard
end


for index=2:size(x,2)-1
    
    
    bla=x; bla(:,[1 index])=[];
    x_ortho = [bla]; 
    res = x(:,index) - x_ortho*pinv(x_ortho)*x(:,index);
    depvar_res_new = y - x_ortho*pinv(x_ortho)*y;
    
    
    
    figure;
    scatter(res, depvar_res_new, 250,'markerfacecolor',[ 0 0 102/256], 'markeredgecolor','k') 
    set(gca,'tickdir', 'out', 'fontname', 'arial')
    xlabel(['Residual Glut/GABA Ratio ',labels_MRS{index-1}], 'fontsize', 95, 'fontname', 'arial')
    ylabel('Residual Effect of DV', 'fontsize',95, 'fontname', 'arial')
    
    ax = gca;
    ax.XColor = 'k';
    ax.YColor = 'k';
    set(gca,'linewidth',2)
    set(gca,'xlim',[-1.9 2.9]);
    
    
    h = lsline;set(h,'color','k','linewidth',3);
    set(gca,'xlim',[-2 3]);
    set(gca,'ylim',[-3 3]);
    L = get(gca,'XLim');
    set(gca,'XTick',[L(1):1:L(2)]);
    
    K = get(gca,'YLim');
    set(gca,'YTick',[K(1):1:K(2)]);
    set(gca,'FontSize',20);
    
    
    [r1,p1,lb,ub]=corrcoef((res), (depvar_res_new));
    
    
    disp('residual correlations:');
    disp([r1(2); p1(2)]); drawnow;
    
    
    if p1(2) <0.07
        keyboard
    end
    
    
    
    % close all
    clear dist res depvar_res_new bla
end

%Resiudalize GABA / glut concentrations for concentrations in all other voxels of interest

for ind=2:(size(x,2))-1
    
    
    check1=x_single; check1(:,[1 ind])=[]; %remove respective GABA concentration from matrix
    x_ortho_gaba = [check1];
    gaba_res = x_single(:,ind) - x_ortho_gaba*pinv(x_ortho_gaba)*x_single(:,ind); %regress those effects out of GABA
    depvar_gaba_res_new = y - x_ortho_gaba*pinv(x_ortho_gaba)*y; %... and out of the dependent variable
    
    clear check1 x_ortho_gaba
    
    index=ind+size(x,2)-2; %do the same for glutamate concentrations
    check1=x_single; check1(:,[1 index])=[];
    x_ortho_glut = [check1];
    glut_res = x_single(:,index) - x_ortho_glut*pinv(x_ortho_glut)*x_single(:,index);
    depvar_glut_res_new = y - x_ortho_glut*pinv(x_ortho_glut)*y;
    
    
    %plot orthogonalized correlations
    figure;
    scatter(gaba_res, depvar_gaba_res_new, 250,'markerfacecolor',[ 0 0 102/256], 'markeredgecolor','k') %; ylim([-1.7 -0.4]); xlim([-0.082 0.105])
    set(gca,'tickdir', 'out', 'fontname', 'arial')
    xlabel(['Residual GABA ',labels_MRS{ind-1}], 'fontsize', 85, 'fontname', 'arial')
    ylabel('Residual Effect of DV', 'fontsize',85, 'fontname', 'arial');
    ax = gca;
    ax.XColor = 'k';
    ax.YColor = 'k';
    set(gca,'linewidth',2)
    
    set(gca,'xlim',[-1.4 1.4]);
    
    h = lsline; set(h,'color','k','linewidth',3);
    set(gca,'xlim',[-1.5 1.5]);
    
    %set actual x and y lim
    set(gca,'xlim',[-3 3]);
    set(gca,'ylim',[-3 3]);
    
    L = get(gca,'XLim');
    set(gca,'XTick',[L(1):1:L(2)])    
    K = get(gca,'YLim');
    set(gca,'YTick',[K(1):1:K(2)]);
    
    
    
    % do the same for glutamate
    figure;
    scatter(glut_res, depvar_glut_res_new, 250,'markerfacecolor',[ 0 0 102/256], 'markeredgecolor','k')%; ylim([-2.4 -0.9]); xlim([-0.2 0.25])
    set(gca,'tickdir', 'out', 'xticklabel', [], 'yticklabel', []);
    xlabel(['Residual Glutamate ',labels_MRS{ind-1}], 'fontsize', 85, 'fontname', 'arial')
    ylabel('Residual Effect of DV', 'fontsize', 95, 'fontname', 'arial')
    
    ax = gca;
    ax.XColor = 'k';
    ax.YColor = 'k';
    set(gca,'linewidth',2)
    ax = gca;
    ax.XColor = 'k';
    ax.YColor = 'k';
    set(gca,'linewidth',2)
    set(gca,'xlim',[-1.4 1.4]);   %set a little shorter xlim to make lsline shorter
    
    h = lsline; set(h,'color','k','linewidth',3);
    
    %set actual x and y lim
    set(gca,'xlim',[-3 3]);
    set(gca,'ylim',[-3 3]);
    
    L = get(gca,'XLim');
    set(gca,'XTick',[L(1):1:L(2)])
    
    K = get(gca,'YLim');
    set(gca,'YTick',[K(1):1:K(2)]);
    
    
    %test their significance
    [r1,p1,lbgaba,ubgaba]=corrcoef((gaba_res), (depvar_gaba_res_new));
    [r2,p2,lbglut,ubglut]=corrcoef((glut_res), (depvar_glut_res_new));
    
    
    disp('residual correlations:');
    disp([r1(2) r2(2); p1(2) p2(2)]); drawnow;
    
    if p1(2) <0.09 || p2(2) <0.09
        keyboard
    end
    
    clear  x_ortho_glut x_ortho_gaba gaba_res glut_res depvar_gaba_res_new depvar_glut_res_new
    close all
end
