
%% simulate and recover artificial data to validate the model
%for questions or criticism,please contact:kaiserl@hhu.de

function[]=sim_model_server_fixtau(iteration,model_idx)

addpath '/gpfs/project/projects/bpsydm/clash/clash_revision/Second_Phase/Modelling/generate_recover/rand_simulate/fixtau'
model={'fixtau_multi'};


num_it=500;

cd('C:\Users\lucak\Desktop\clash_final\logfiles\');

vp = {'1'; '2'; '3'; '4'; '5'; '7'; '9'; '10'; '12';  '14'; '15'; '16'; '17'; '18';  '20';  '22'; '23'; '24';  '26'; '27'; '28'; '29'; '30'; '31'; '32'; '33'; '34'; '35'; '38';};

tau=6.62;

disp(model{model_idx});
disp(vp{iteration});
data = importdata(['clashMEG',vp{iteration},'.txt']);

data = data.data;

data(data(:,5)== 999,:)=[];        % no breaks
org_data = data;


% call model estimation
[distances] = main_model_SU_ourdistort_sim_fixtau(org_data,model{model_idx},'lu',num_it,vp{iteration},tau);
cp_fixtau.distances = distances;



cd    (['/gpfs/project/projects/bpsydm/clash/clash_revision/Second_Phase/Modelling/generate_recover/' model{model_idx} 'SU/'])


save(['sim_ourdistort_fixtau' model{model_idx}  '_' num2str(iteration) 'lu'], 'cp_fixtau');

end

