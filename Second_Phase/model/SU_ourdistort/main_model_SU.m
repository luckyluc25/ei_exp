function [fitparams,betas, modelfit] = main_model_SU(data,model,dat_class,num_it,tau)

lb =    [0 0 0 0 0 0 0 ];     % lower bound constraint

for start=1:num_it
    
    if strcmp(model,'multi')
        ub=[60 3 3];
        x0 = [ub(1)*rand(1) 3*rand(1) 3*rand(1)];
    elseif strcmp(model,'adi')
        ub=[60 60 3 3 ];
        x0 = [ub(1)*rand(1) ub(2)*rand(1) 3*rand(1) 3*rand(1)];
    elseif strcmp(model,'hybrid') || strcmp(model,'aditau')
        ub=[60 60 60 3 3];
        x0 = [ub(1)*rand(1) ub(2)*rand(1) ub(3)*rand(1) 3*rand(1) 3*rand(1)];
    elseif strcmp(model,'fixtau_multi')
        ub=[3 3];
        x0 = [3*rand(1) 3*rand(1)];
    end
    
    options = optimoptions('fmincon','Display','off');
    [fitparams_between, modelfit_between] = fmincon(@model_distort, x0, [], [], [], [], lb,ub,[],options);
    fit_all(start)=modelfit_between;
    fitparams_all(start,:)=fitparams_between;
    
end

minimum = min(fit_all);
[x]=find(fit_all==minimum);

if length(x)>1
    rand_idx=randi(length(x));
    x=x(rand_idx);
end

modelfit=fit_all(x);
fitparams=fitparams_all(x,:);

if strcmp(model,'hybrid') || strcmp(model, 'hybridtau')
    beta_mult_all=(fitparams(1)/(fitparams(1)+fitparams(2)+fitparams(3)));
    beta_m_all=fitparams(2)/(fitparams(2)+fitparams(3));
    betap_all=fitparams(3)/(fitparams(2)+fitparams(3));
    betas=[beta_mult_all beta_m_all betap_all];
    
elseif strcmp(model,'adi')|| strcmp(model,'aditau')
    
    beta_m_all=fitparams(1)/(fitparams(1)+fitparams(2));
    betap_all=fitparams(2)/(fitparams(1)+fitparams(2));
    
    betas=[beta_m_all betap_all];
    
else
    betas=[];
    
end

%=============== all done ========================= %


% ================ PT model ========================= %
    function energy = model_distort(params)
        
        
        if strcmp(dat_class,'lu')
            
            choices=data(:,13);
            
            rp_l = data(:,10)/10;
            rp_r = data(:,12)/10;
            
            %rescale mags to make sure they are in the same range for every
            %participant
            rm_all=[data(:,9) data(:,11)];
            rm_all=rescale(rm_all,1,10);
            rm_l=rm_all(:,1);
            rm_r=rm_all(:,2);
            
            
        elseif strcmp(dat_class,'ger')
            choices=data(:,11);
            rp_l = data(:,3)/100;
            rp_r = data(:,4)/100;
            
            
            rm_all=[data(:,5) data(:,6)];
            rm_all=rescale(rm_all,1,10);
            rm_l=rm_all(:,1);
            rm_r=rm_all(:,2);
            
            
        end
        
        alpha=params(end-1);
        gamma=params(end);
        
        sm_l = rm_l.^alpha;
        sm_r = rm_r.^alpha;
        sp_l = (rp_l .^gamma) ./ ( (rp_l.^gamma  +  (1-rp_l).^gamma).^(1./gamma)  );
        sp_r = (rp_r .^gamma) ./ ( (rp_r.^gamma  +  (1-rp_r).^gamma).^(1./gamma)  );
        
        
        
        
        clear rm_r rm_l rp_l rp_r
        if strcmp(model,'multi')
            omegamult = params(1);
            
            SV_l=(omegamult)*(sp_l.*sm_l);
            SV_r=(omegamult)*(sp_r.*sm_r);
            
            
        elseif strcmp(model,'fixtau_multi')
            omegamult = tau;
            
            SV_l=(omegamult)*(sp_l.*sm_l);
            SV_r=(omegamult)*(sp_r.*sm_r);
            
        elseif strcmp(model,'adi')
            omegam = params(1);
            omegap=params(2);
            
            SV_l=(omegam*sm_l+omegap*sp_l);
            SV_r=(omegam*sm_r+omegap*sp_r);
            
        elseif strcmp(model,'hybrid')
            omegamult = params(1);
            omegam = params(2);
            omegap=params(3);
            
            beta_mult=round((omegamult/(omegam+omegap+omegamult)),4);
            betam=round(omegam/(omegam+omegap),4);
            betap=round(omegap/(omegam+omegap),4);
            
            SV_l=((omegam+omegap+omegamult))*((1-beta_mult)*(betam*sm_l+betap*sp_l)+beta_mult*(sp_l.*sm_l));
            SV_r=((omegam+omegap+omegamult))*((1-beta_mult)*(betam*sm_r+betap*sp_r)+beta_mult*(sp_r.*sm_r));
            
            
        end
        
        
        
        
        
        pch(:,1) =  1./(1 + exp(-(SV_l-SV_r)));     % modeled probability for left choice
        pch(:,2) =1./(1 + exp(-(SV_r-SV_l))) ;
        
        subpch=zeros(length(choices),1);
        for i=1:length(pch)
            subpch(i,1) = pch(i,choices(i));
        end
        subpch(subpch<10^-100) = 10^-100;
        energy = -sum(log(subpch));
        
        
    end
% ==================================================== %


end     % end of parent function









