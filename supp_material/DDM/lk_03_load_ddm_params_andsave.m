clear; clc;close all
model_int={'standard_ddm_var_bias_nosv'};

for modix=1:numel(model_int)
    
    
    cd (['C:\Users\lucak\Desktop\clash_script_upload\DDM_FirstPhase\short_results\' model_int{modix}]);
    gel_rub=readtable('gelman_rubin.txt');
    
    diff=sum(abs(1-gel_rub.Var2)>0.02);
    
    %check whether chains have not converged
    if diff~=0
        keyboard
        index=find(abs(1-gel_rub.Var2)>0.02);
        index_test=gel_rub.Var1(index);
    end
    
    %check results
    
    check=readtable(['results-combined.csv']);
    search=check.Var1;
    
    %set up empty variables
    a_single=[]
    v_single=[]
    t_single=[]
    z_single=[]
    
    a_intercept=[];
    t_intercept=[];
    v_intercept=[];
    z_intercept=[];
    
    
    %a
    
    index = find(contains(search,'a_subj'));
    if ~isempty(index)
        a_single=check(index(1:29),2);
        a_single=a_single.mean;
    end
    
    clear index
    
    index = find(contains(search,'a_Intercept_su'));
    a_intercept=check(index,2);
    a_intercept=a_intercept.mean;
    
    clear index
    
    index = find(contains(search,'a_costs_su'));
    a_costs=check(index,2);
    a_costs=a_costs.mean;
    
    clear index
    
    
    
    %z
    
    index = find(contains(search,'z_Intercept_subj'));
    z_intercept=check(index,2)
    z_intercept=z_intercept.mean;
    
    clear index
    
    index = find(contains(search,'z_subj'));
    if ~isempty(index) && isempty(z_intercept)
        z_single=check(index,2);
        z_single=z_single.mean;
    end
    
    clear index
    
    
    %t
    clear index
    
    index = find(contains(search,'t_Intercept_su'));
    if ~isempty(index)
        t_intercept=check(index,2);
        t_intercept=t_intercept.mean;
    end
    
    clear index
    
    index = find(contains(search,'t_subj'));
    if ~isempty(index) && isempty(t_intercept)
        t_single=check(index(1:29),2);
        t_single=t_single.mean;
    end
    
    clear index
    
    
    
    
    %v
    index = find(contains(search,'v_subj'));
    if ~isempty(index)
        v_single=check(index,2);
        v_single=v_single.mean;
    end
    
    clear index
    
    index = find(contains(search,'v_Intercept_su'));
    v_intercept=check(index,2);
    v_intercept=v_intercept.mean;
    
    clear index
    
    
    
    index = find(contains(search,'v_costs_su'));
    v_costs=check(index,2);
    v_costs=v_costs.mean;
    
    clear index
    
    index = find(contains(search,'v_val_sum_subj'));
    v_valsum=check(index,2);
    v_valsum=v_valsum.mean;
    
    clear index
    
    
    
    model.a_single=a_single;
    model.t_single=t_single;
    model.v_single=v_single;
    model.z_single=z_single;
    
    model.a_intercept=a_intercept;
    model.v_intercept=v_intercept;
    model.t_intercept=t_intercept;
    model.z_intercept=z_intercept;
    
    
    
    model.a_costs=a_costs;
    model.v_costs=v_costs;
    
    
    
    save(['C:\Users\lucak\Desktop\clash_script_upload\First_Phase\DV\' model_int{modix}],'model')
    
    %plot posteriors
    nodes=importdata('group_traces.csv');
    
    for nodes_idx=2:length(nodes.textdata)
        figure()
        hist(nodes.data(:,nodes_idx));
        title(nodes.textdata(nodes_idx));
        figure()
        plot(nodes.data(:,nodes_idx));
        title(nodes.textdata(nodes_idx));
        
    end
    
       
end




