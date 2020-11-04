cd('C:\Users\lucak\Desktop\clash_script_upload\Second_Phase\DV')
per_correct=importdata('objcorP.mat');
rt_2=importdata('rt_2nd.mat');
rt_effects=importdata('c_feb.mat');
val_diff_rt_vg=rt_effects(:,2);


%modeling effects
win_model=importdata('C:\Users\lucak\Desktop\clash_final\Second_Phase\Modelling_lokal\extended_fitting\SU_ourdistortfixtau_multilu.mat');


corr_mat=[per_correct rt_2 val_diff_rt_vg win_model.params(:,1) win_model.params(:,2)];

check=corr(corr_mat);
imagesc(check);


