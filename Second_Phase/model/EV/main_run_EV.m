


%% run all EV models
%for questions or criticism,please contact:kaiserl@hhu.de

addpath  C:\Users\lucak\Desktop\clash_final\Second_Phase\Modelling_lokal
addpath C:\Users\lucak\Desktop\clash_final\Second_Phase\Modelling_lokal\EV
cd('C:\Users\lucak\Desktop\clash_final\logfiles\');

model={'multi','adi','hybrid'}; %,'fixtau_multi'};
vp = {'1'; '2'; '3'; '4'; '5'; '7'; '9'; '10'; '12';  '14'; '15'; '16'; '17'; '18';  '20';  '22'; '23'; '24';  '26'; '27'; '28'; '29'; '30'; '31'; '32'; '33'; '34'; '35'; '38';};
num_it=1000;

for model_idx=1:numel(model)
    cp_linear.names = {'omegam', 'omegap', 'omegamult','LLE', 'bic' };
    vp = {'1'; '2'; '3'; '4'; '5'; '7'; '9'; '10'; '12';  '14'; '15'; '16'; '17'; '18';  '20';  '22'; '23'; '24';  '26'; '27'; '28'; '29'; '30'; '31'; '32'; '33'; '34'; '35'; '38';};
    cd('C:\Users\lucak\Desktop\clash_script_upload\logfiles\');
    
    for iteration=1:length(vp)
        disp(model{model_idx});
        disp(vp{iteration});
        data = importdata(['clashMEG',vp{iteration},'.txt']);
        
        data = data.data;
        
        data(data(:,5)== 999,:)=[];        % no breaks
        org_data = data;
        
        if strcmp(model{model_idx},'fixtau_multi')
            parameter=importdata('C:\Users\lucak\Desktop\clash_final\Second_Phase\Modelling_lokal\extended_fitting\EVmultilu.mat');
            parameter=parameter.params;
            med=median(parameter);
            tau=med(1);
        else
            tau=[];
        end
        tau=round(tau,2);
        
        % call model estimation
        [params, betas,fitval] = main_model_EV(org_data,model{model_idx},'lu',num_it,tau);
        bic = 2*fitval + numel(params)*log(size(data,1));
        aic=2*fitval + (2*numel(params));
        cp_EV.params(iteration,:) = [params, betas,fitval, bic,aic];
        
    end
    
    cd 'C:\Users\lucak\Desktop\clash_final\Second_Phase\Modelling_lokal\extended_fitting'
    
    save(['EV' model{model_idx} 'lu'], 'cp_EV');
    
    clear cp_EV bic aic params fitval betas
    
    cp_linear.names = {'omegam', 'omegap', 'omegamult','LLE', 'bic' };
    vp = {'01'; '02'; '03'; '04'; '05'; '06'; '07'; '08'; '09';  '10'; '11'; '12'; '13'; '14';  '15';  '16'; '17'; '18';  '19'; '20'; '21'; '22'; '23'; '24'; '25'};
    cd('C:\Users\lucak\Desktop\clash_revision\gerhard_log\');
    %
    %
    for iteration=1:length(vp)
        disp(vp{iteration});
        data = importdata(['vp',vp{iteration},'_freechoice.txt']);
        
        
        org_data = data;
        
        
        if strcmp(model{model_idx},'fixtau_multi')
            parameter=importdata('C:\Users\lucak\Desktop\clash_final\Second_Phase\Modelling_lokal\extended_fitting\EVmultiger.mat');
            parameter=parameter.params;
            med=median(parameter);
            tau=med(1);
        else
            tau=[];
        end
        tau=round(tau,2);
        
        % call model estimation
        [fitparams,betas, modelfit] = main_model_EV(data,model{model_idx},'ger',num_it);
        bic = 2*modelfit + numel(fitparams)*log(size(data,1));
        aic=2*modelfit + (2*numel(fitparams));
        cp_EV.params(iteration,:) = [fitparams, betas,modelfit,bic,aic];
        
    end
    
    cd   'C:\Users\lucak\Desktop\clash_final\Second_Phase\Modelling_lokal\extended_fitting'
    save(['EV' model{model_idx} 'ger'], 'cp_EV');
    
    clear cp_EV
end

