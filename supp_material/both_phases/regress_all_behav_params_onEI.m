%clean up
clear; close all; clc 


%add the paths you need


cd( 'C:\Users\lucak\Desktop\clash_script_upload\First_Phase\DV');

files=dir('*.mat')
for k=1:numel(files)
    load(files(k).name);
end

%load dependent variables
cd( 'C:\Users\lucak\Desktop\clash_script_upload\Second_Phase\DV');

files=dir('*.mat')
for k=1:numel(files)
    load(files(k).name);
end


%% read in GABA /glutamate concentrations
vp = {'1'; '2'; '3'; '4'; '5'; '7'; '9'; '10'; '12';  '14'; '15'; '16'; '17'; '18';  '20';  '22'; '23'; '24';  '26'; '27'; '28'; '29'; '30'; '31'; '32'; '33'; '34'; '35'; '38'};

table=importdata( 'C:\Users\lucak\Desktop\clash_script_upload\ressources_share\table_anonym.mat');

%rename combines MR names with subject numbers
rename=importdata( 'C:\Users\lucak\Desktop\clash_script_upload\ressources_share\rename_anonym.mat');

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
        
        %spm12 idx is the index of the subject; 
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

x = [normalise(x)];
model= importdata('C:\Users\lucak\Desktop\clash_final\Second_Phase\Modelling_lokal\extended_fitting\SU_ourdistortfixtau_multilu.mat');
design=[switch_index linReg_1stphase_check(:,4) (rt_1st)' (rt_2nd) objcorP c_feb(:,2) model.params(:,1) model.params(:,2)];
design=normalise(design);


vmPFC=x(:,4);
dACC=x(:,5);    


lm = fitlm(design,vmPFC);
ls = fitlm(design,dACC);
coefCI(ls)

