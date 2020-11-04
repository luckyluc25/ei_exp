
function [distances] = main_model_SU_ourdistort_sim_fixtau(data,model,dat_class,num_it,vp,tau)

lb =    [0 0 0 0 0 0 0 ];     % lower bound constraint
ub =    [60 60 60 60 60 60];   % upper bound constraint

for start=1:num_it
    
    if strcmp(model,'multi')
        ub=[60 3 3];
        x0 = [ub(1)*rand(1) 3*rand(1) 3*rand(1)];
    elseif strcmp(model,'adi')
        ub=[60 60 3 3];
        x0 = [ub(1)*rand(1) ub(2)*rand(1) ub(3)*rand(1) ub(4)*rand(1)];
    elseif strcmp(model,'hybrid') || strcmp(model,'aditau')
        ub=[60 60 60 3 3];
        x0 = [ub(1)*rand(1) ub(2)*rand(1) ub(3)*rand(1) ub(4)*rand(1) ub(5)*rand(1)];
    elseif strcmp(model,'hybridtau')
        x0 = [ub(1)*rand(1) ub(2)*rand(1) ub(3)*rand(1) ub(4)*rand(1) ub(5)*rand(1) ub(6)*rand(1)];
    elseif strcmp(model,'fixtau_multi')
        ub=[ 3 3];
        x0 = [3*rand(1) 3*rand(1)];
    end
    
    options = optimoptions('fmincon','Display','off');
       
    %start simulation here
    
    [rev_params,rev_betas]=make_data(data,model,x0,num_it,dat_class,lb,ub,tau);
    
    distances(start,:)=[x0,rev_params,rev_betas];   
           
end



%=============== all done ========================= %



end     % end of parent function









