
function [fitparams,betas, modelfit] = main_model_EV(data,model,dat_class,num_it,tau)

lb =    [0 0 0 0 0 0];     % lower bound constraint
ub =    [60 60 60 60 60];   % upper bound constraint


for start=1:num_it
    
    if strcmp(model,'multi')
        x0 = [ub(1)*rand(1)];
    elseif strcmp(model,'adi')
        x0 = [ub(1)*rand(1) ub(2)*rand(1)];
    elseif strcmp(model,'hybrid') || strcmp(model,'aditau')
        x0 = [ub(1)*rand(1) ub(2)*rand(1) ub(3)*rand(1)];
    elseif strcmp(model,'hybridtau')
        x0 = [ub(1)*rand(1) ub(2)*rand(1) ub(3)*rand(1) ub(4)*rand(1)];
    end
    
    options = optimoptions('fmincon','Display','off');
    [fitparams_between, modelfit_between] = fmincon(@lin_mod, x0, [], [], [], [], lb,ub,[],options);
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
    
    betas=[ beta_m_all betap_all];
    
else
    betas=[];
    
end

%=============== all done ========================= %


% ================ PT model ========================= %
    function energy = lin_mod(params)
        
        
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
        
        
        if strcmp(model,'multi')
            omegamult = params(1);
            
            SV_l=(omegamult)*(rp_l.*rm_l);
            SV_r=(omegamult)*(rp_r.*rm_r);
            
        elseif strcmp(model,'adi')
            omegam = params(1);
            omegap=params(2);
            
            SV_l=(omegam*rm_l+omegap*rp_l);
            SV_r=(omegam*rm_r+omegap*rp_r);
            
        elseif strcmp(model,'hybrid')
            omegamult = params(1);
            omegam = params(2);
            omegap=params(3);
            
            beta_mult=round((omegamult/(omegam+omegap+omegamult)),4);
            betam=round(omegam/(omegam+omegap),4);
            betap=round(omegap/(omegam+omegap),4);
            
            SV_l=((omegam+omegap+omegamult))*((1-beta_mult)*(betam*rm_l+betap*rp_l)+beta_mult*(rp_l.*rm_l));
            SV_r=((omegam+omegap+omegamult))*((1-beta_mult)*(betam*rm_r+betap*rp_r)+beta_mult*(rp_r.*rm_r));
            
            
        elseif strcmp(model,'fixtau_multi')
            omegamult = tau;
            
            SV_l=(omegamult)*(rp_l.*rm_l);
            SV_r=(omegamult)*(rp_r.*rm_r);
            
        elseif strcmp(model,'aditau')
            omegam = params(1);
            omegap=params(2);
            tau=params(3);
            
            SV_l=(tau)*(round(omegam/(omegam+omegap),4)*rm_l+round(omegap/(omegam+omegap),4)*rp_l);
            SV_r=(tau)*(round(omegam/(omegam+omegap),4)*rm_r+round(omegap/(omegam+omegap),4)*rp_r);
            
        elseif strcmp(model, 'hybridtau')
            omegamult = params(1);
            omegam = params(2);
            omegap=params(3);
            tau=params(4);
            
            
            beta_mult=round((omegamult/(omegam+omegap+omegamult)),4);
            betam=round(omegam/(omegam+omegap),4);
            betap=round(omegap/(omegam+omegap),4);
            
            SV_l=((tau))*((1-beta_mult)*(betam*rm_l+betap*rp_l)+beta_mult*(rp_l.*rm_l));
            SV_r=((tau))*((1-beta_mult)*(betam*rm_r+betap*rp_r)+beta_mult*(rp_r.*rm_r));
        end
        
        
        
        
        
        pch(:,1) =  1./(1 + exp(-(SV_l-SV_r)));
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









