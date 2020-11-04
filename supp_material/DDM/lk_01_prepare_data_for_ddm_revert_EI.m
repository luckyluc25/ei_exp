%% make a table for DDM model fitting

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
% 10=   2nd stage: reward probability left
% 11=   2nd stage: reward magnitude right 
% 12=   2nd stage: reward probability right 
% 13=   choice (left/right) at 2nd stage
% 14=   RT for 2nd stage choice
% 15=   outcome (reward/no reward) left for 2nd stage
% 16=   outcome (reward/no reward) right for 2nd stage
% 17=   cumulative points won by subject

%=========================================================================


%% read in GABA /glutamate concentrations
vp = {'1'; '2'; '3'; '4'; '5'; '7'; '9'; '10'; '12';  '14'; '15'; '16'; '17'; '18';  '20';  '22'; '23'; '24';  '26'; '27'; '28'; '29'; '30'; '31'; '32'; '33'; '34'; '35'; '38'};

% E/I balances are read out from csv files to a table and stored under % their MR name
table=importdata(  'C:\Users\lucak\Desktop\clash_script_upload\ressources_share\table_anonym.mat');

%rename combines MR names with subject numbers
rename=importdata( 'C:\Users\lucak\Desktop\clash_script_upload\ressources_share\rename_anonym.mat');

%in gm_spm12_test grey matter ceoncentrations are stored by summing across %P(GM). They are later divided by total number of voxels in mask to
%approximate relative GM concentrations.
gm_spm12_gmwm=importdata('C:\Users\lucak\Desktop\clash_script_upload\ressources_share\gm_spm12_gmwm_anonym.mat') ;

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

valdiff_all=[];
pvd_all=[];
matrix_all=[];

for iteration=1:numel(vp)
    
    data=importdata(['C:\Users\lucak\Desktop\clash_script_upload\logfiles\clashMEG',vp{iteration},'.txt']);
    data=data.data;
    data=data(data(:,5)~=999,:);   
    ntrl=size(data,1);
    
    rp=data(:,[10 12])/10;
    rm=data(:,[9 11]);
    
    ch_side=data(:,5); ch_side=[0; abs(diff(ch_side))];
    
    rt=data(:,8);
    sw = [0; abs(diff(data(:,6)))];
    
    
    sw(sw==1)=121;
    sw(sw==0)=1;
    sw(sw==121)=0;
    
    wins=zeros(size(rt));
    wins( (data(:,13)==1 & data(:,15)==1) | (data(:,13)==2 & data(:,16)==1) ) = 1;       % did subject win at 2nd stage?
    
    k=3;                   
    warning('off');
    repwins=toeplitz(wins,zeros(1,k+1));
    warning('on');
    
    repwins=repwins(:,2:end);     
    
    
    vcurr=zeros(ntrl,1);                    
    valt=zeros(ntrl,1);
    
    for k=2:ntrl                           
        if data(k-1,6)==1                                           
            vcurr(k) = data(k-1,3);
            valt(k) = data(k-1,4);
        elseif data(k-1,6)==2
            vcurr(k) = data(k-1,4);
            valt(k) = data(k-1,3);
        end
    end
    
    pvd=valt-vcurr;
    costs=data(:,7);
    trial_number=1:length(costs);
    
    out_idx=find(rt>0.3& (rt<4));
    trl_idx=out_idx;
   
    my_matrix= [repmat(iteration,length(pvd(trl_idx)),1) rt(trl_idx)  sw(trl_idx) zscore(pvd(trl_idx))  (costs(trl_idx)) (trial_number(trl_idx)') ch_side(trl_idx) repwins(trl_idx)];
    
    EI_subj=EI(iteration,:);EI_subj=repmat(EI_subj,length(my_matrix),1);
    my_matrix=[my_matrix EI_subj];
    
    matrix_all=[matrix_all;my_matrix] ;
    
    
    
   
    
end

matrix_all(:,[5 6])=zscore(matrix_all(:,[5 6]));

T=array2table(matrix_all);
T.Properties.VariableNames(1:13) = {'subj_idx','rt','response','pvd','costs','trial_number','ch_side','repwins','dlpfc','m1li','m1re','vmpfc','dacc'};

cd ( 'C:\Users\lucak\Desktop\');
writetable(T,'hddm_clash_patch_revert_EI.csv','Delimiter','comma');


