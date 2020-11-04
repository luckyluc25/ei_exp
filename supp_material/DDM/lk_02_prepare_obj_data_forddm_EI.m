%% estimate the effects of decision variables on log RT

%for questions or criticism, contact kaiserl@hhu.de

%%02/20 added choice vector to design matrix

%=========================================================================
% required input: behavioral data (.txt-file with column headers), columns have are organized:
% 1 =   subject ID
% 2 =   trial no
% 3 =   value of patch A
% 4 =   value of patch B
% 5 =   position of patches on screen (1=A left, B right, 2=A right, B left) Caution: 999=Pause   
% 6 =   patch choice on current trial (
% 7 =   switch cost
% 8 =   RT for patch decision
% 9 =   2nd stage: reward magnitude left
% 10=  2nd stage: reward probability left
% 11=  2nd stage: reward magnitude right 
% 12=  2nd stage: reward probability right 
% 13=  choice (left/right) at 2nd stage
% 14=  RT for 2nd stage choice
% 15=  outcome (reward/no reward) left for 2nd stage
% 16=  outcome (reward/no reward) right for 2nd stage
% 17=  cumulative points won by subject

%=========================================================================
clear all; close all;

%%falls wir etwas ändern wäre es vermutlich besser EI am Ende zu scoren?!

vp = {'1'; '2'; '3'; '4'; '5'; '7'; '9'; '10'; '12';  '14'; '15'; '16'; '17'; '18';  '20';  '22'; '23'; '24';  '26'; '27'; '28'; '29'; '30'; '31'; '32'; '33'; '34'; '35'; '38'};

% E/I balances are read out from csv files to a table and stored under % their MR name
table=importdata( 'C:\Users\lucak\Desktop\clash_share-main\ressources_share\table_anonym.mat');

%rename combines MR names with subject numbers
rename=importdata( 'C:\Users\lucak\Desktop\clash_share-main\ressources_share\rename_anonym.mat');

%in gm_spm12_test grey matter ceoncentrations are stored by summing across %P(GM). They are later divided by total number of voxels in mask to
%approximate relative GM concentrations.
gm_spm12_gmwm=importdata('C:\Users\lucak\Desktop\clash_script_upload\ressources_share\gm_spm12_gmwm_anonym.mat') ;


%just a check-if one variable is not there that should be, the code will
%break!
gm_spm12_gmwm.gm_all(gm_spm12_gmwm.gm_all==0)=NaN;
gm_spm12_gmwm.wm_all(gm_spm12_gmwm.wm_all==0)=NaN;

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
        
       
        %where is the correct GM concenctration stored?
        spm12_idx=(~cellfun(@isempty,(strfind(gm_spm12_gmwm.mrName,mrName))));
        concen_con=find((~cellfun(@isempty,(strfind(erase(gm_spm12_gmwm.vox_name(spm12_idx,:),'_'),name))))); 
        
        
        if isempty (idx)
            glut_corr(count,concen_con)=NaN;
            gaba_corr(count,concen_con)=NaN;
            count=count+1;
            continue
        end
        
        %Normalise by GM within each voxel    
        
        %spm12 idx is the index of the subject; concentrations are stored in the same order as they are read out here
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


x=[(glut_corr./gaba_corr)];
EI=normalise(x);

%add the path you need
addpath('C:\Users\luca\Desktop\magma\MatFunc');
valdiff_all=[];
pvd_all=[];
matrix_all=[];

vp = {'1'; '2'; '3'; '4'; '5'; '7'; '9'; '10'; '12';  '14'; '15'; '16'; '17'; '18';  '20';  '22'; '23'; '24';  '26'; '27'; '28'; '29'; '30'; '31'; '32'; '33'; '34'; '35'; '38';};

for iteration=1:numel(vp)
    
    % load in data
    data=importdata(['C:\Users\lucak\Desktop\clash_share-main\logfiles\clashMEG',vp{iteration},'.txt']);
    data=data.data;
    data=data(data(:,5)~=999,:);     % get rid of pauses
    
    %rescale them as they have been rescaled in modelling
    rp=data(:,[10 12])/10;
    rm=data(:,[9 11]);
    
    
    ev=rp.*rm;
    %load distortion params
    
    choices=data(:,13);
    rt=data(:,14);
    sw = [0; abs(diff(data(:,6)))];
    
   ntrl=size(data,1);
    
    for k=1:ntrl
        if data(k,6)==1
            
            vcurr(k) = data(k,3);
            valt(k) = data(k,4);
        elseif data(k,6)==2
            vcurr(k) = data(k,4);
            valt(k) = data(k,3);
        end
    end
    
    pvd=valt-vcurr;
    
    
    nb_obj = double((rp(:,1)>=rp(:,2) & rm(:,1)>=rm(:,2)) | (rp(:,2)>=rp(:,1) & rm(:,2)>=rm(:,1)));
    nb_obj(nb_obj==1)=2;
    nb_obj(nb_obj==0)=1;
    
    
    trial_number=1:length(nb_obj);
    
       
    %code responses according to accuracy
    all=find(ev(:,1)~=ev(:,2));    
    maxev=double(ev(:,1)>ev(:,2));
    maxev(ev(:,1)<ev(:,2))=2;
    
    objcorr=(choices==maxev);
    
    out_idx=find(rt<0.3| (rt>4));
    
    if ~isempty(out_idx)
        all(ismember(all,out_idx))=[];
    end
    
    val_diff=abs(ev(:,1)-ev(:,2)); 
    val_sum=abs(ev(:,1)+ev(:,2));
    
    
    costs=data(:,7);
    
    index=3;
    choice_right=double(choices==1);
    repchoices= toeplitz(choice_right',zeros(1,index+1));
    repchoices=repchoices(:,2);
    
    all_wins=double((data(:,16)==1 & data(:,13)==2) | (data(:,15) & data(:,13)==1));
    all_wins=toeplitz(all_wins,zeros(1,index+1));
    all_wins=all_wins(:,2);
    
    choices_pvd=zeros(length(choices),1);
    
    for k_idx=1:length(choices)
        if data(k_idx,5)==1
            if data(k_idx,6)==1
                choices_pvd(k_idx)=1;
            elseif data(k_idx,6)==2
                choices_pvd(k_idx)=2;
            end
        elseif data(k_idx,5)==2
            if data(k_idx,6)==1
                choices_pvd(k_idx)=2;
            elseif data(k_idx,6)==2
                choices_pvd(k_idx)=1;
            end
        end
    end
    
    
    my_matrix= [repmat(iteration,length(all),1) rt(all)  objcorr(all) zscore(val_diff(all)) zscore(val_sum(all)) nb_obj(all) zscore(pvd(all))'  (trial_number(all)') sw(all) (costs(all)) repchoices(all) all_wins(all) choices_pvd(all) choices(all)];
    
    EI_subj=EI(iteration,:);EI_subj=repmat(EI_subj,length(my_matrix),1);
    my_matrix=[my_matrix EI_subj];
    
    
    matrix_all=[matrix_all;my_matrix] ;
    
    
    
    
    
end

matrix_all(:,[8 10])=zscore(matrix_all(:,[ 8 10]));

T=array2table(matrix_all);
T.Properties.VariableNames(1:19) = {'subj_idx','rt','response','val_diff','val_sum','nb','pvd','trial_number','sw','costs','repchoices','all_wins','choices_pvd','choices','dlpfc','m1li','m1re','vmpfc','dacc'};

cd ( 'C:\Users\lucak\Desktop\');
writetable(T,'hddm_obj_clash_EI.csv','Delimiter','comma');



