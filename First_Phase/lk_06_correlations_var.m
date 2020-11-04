
cd('C:\Users\lucak\Desktop\clash_script_upload\First_Phase\DV');

pla=importdata('switch_index.mat');
rt_effects=importdata('linReg_1stphase_check.mat');
cost_rt=rt_effects(:,4);
rt_1=importdata('rt_1st.mat');


corr_mat=[pla (rt_1)' cost_rt];
check=corr(corr_mat);
imagesc(check);


