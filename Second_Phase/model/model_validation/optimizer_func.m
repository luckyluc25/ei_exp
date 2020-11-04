function energy = optimizer_func(params,data,dat_class,model,tau)


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
    
    
elseif strcmp(model,'fixtau_multi')
    omegamult = tau;
    
    SV_l=(omegamult)*(sp_l.*sm_l);
    SV_r=(omegamult)*(sp_r.*sm_r);
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

