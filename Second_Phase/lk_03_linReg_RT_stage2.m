%% estimate the effects of decision variables on log RT

%for questions or criticism, contact kaiserl@hhu.de

%=========================================================================
% required input: behavioral data (.txt-file with column headers), columns have are organized:
% 1 =   subject ID
% 2 =   trial no
% 3 =   value of patch A
% 4 =   value of patch B
% 5 =   position of patches on screen (1=A left, B right, 2=A right, B left) Caution: 999=Pause
% 6 =   patch choice on current trial (
% 7 =   switch cost
% 8 =   RT for patch decision
% 9 =   2nd stage: reward magnitude left
% 10=  2nd stage: reward probability left
% 11=  2nd stage: reward magnitude right
% 12=  2nd stage: reward probability right
% 13=  choice (left/right) at 2nd stage
% 14=  RT for 2nd stage choice
% 15=  outcome (reward/no reward) left for 2nd stage
% 16=  outcome (reward/no reward) right for 2nd stage
% 17=  cumulative points won by subject

%=========================================================================
clear all; close all;

%add the path you need
addpath('C:\Users\luca\Desktop\magma\MatFunc');
addpath('C:\Users\lucak\Desktop\measures-of-effect-size-toolbox-master');


vp = {'1'; '2'; '3'; '4'; '5'; '7'; '9'; '10'; '12';  '14'; '15'; '16'; '17'; '18';  '20';  '22'; '23'; '24';  '26'; '27'; '28'; '29'; '30'; '31'; '32'; '33'; '34'; '35'; '38';};

for iteration=1:numel(vp)
    
    data=importdata(['C:\Users\lucak\Desktop\clash_script_upload\logfiles\clashMEG',vp{iteration},'.txt']);
    data=data.data;
    data=data(data(:,5)~=999,:);     % get rid of pauses
    
    rp=data(:,[10 12])/10;
    rm=data(:,[9 11]);
    
    choices=data(:,13);
    rt=data(:,14);
    rt_2nd(iteration,:)=median(rt);
    sw = [0; abs(diff(data(:,6)))];
    
    ntrl=size(data,1);
    
    %We use the patch value differences of the current trial here, since the
    %participant already knows them during value guided choice
    for k=1:ntrl
        if data(k,6)==1
            vcurr(k) = data(k,3);
            valt(k) = data(k,4);
        elseif data(k,6)==2
            vcurr(k) = data(k,4);
            valt(k) = data(k,3);
        end
    end
    
    
    pvd=valt-vcurr; costs=data(:,7);
    
    
    ev=rp.*rm;
    
    choices_pvd=zeros(length(choices),1);
    for k_idx=1:length(choices)
        if data(k_idx,5)==1
            if data(k_idx,6)==1
                choices_pvd(k_idx)=1;
            elseif data(k_idx,6)==2
                choices_pvd(k_idx)=2;
            end
        elseif data(k_idx,5)==2
            if data(k_idx,6)==1
                choices_pvd(k_idx)=2;
            elseif data(k_idx,6)==2
                choices_pvd(k_idx)=1;
            end
        end
    end
    
    choices_pvd(choices_pvd==1) = 0; choices_pvd(choices_pvd==2) = 1;
    
    index=3;
    choice_right=double(choices==2);
    repchoices= toeplitz(choice_right',zeros(1,index+1));
    repchoices=repchoices(:,2:end);
    
    
    all_wins=double((data(:,16)==1 & data(:,13)==2) | (data(:,15) & data(:,13)==1));
    all_wins=toeplitz(all_wins,zeros(1,index+1));
    all_wins=all_wins(:,2:end);
    
    
    nb_obj = double((rp(:,1)>=rp(:,2) & rm(:,1)>=rm(:,2)) | (rp(:,2)>=rp(:,1) & rm(:,2)>=rm(:,1)));
    
    
    sum_nb_obj(iteration,:)=sum(nb_obj)/length(nb_obj);
    
    
    trial_count=[1:length(nb_obj)]';
    
    x=zscore([(abs(ev(:,1)-ev(:,2))) ev(:,1)+ev(:,2) costs (pvd') repchoices(:,1) all_wins(:,1) choices_pvd nb_obj trial_count sw]);
    
    [linReg2]=fitlm(x,zscore(log(rt)));
    
    c_feb(iteration,:)=linReg2.Coefficients.Estimate;
    
    
end

save( 'C:\Users\lucak\Desktop\clash_script_upload\Second_Phase\DV\rt_2nd','rt_2nd');
save( 'C:\Users\lucak\Desktop\clash_script_upload\Second_Phase\DV\c_feb','c_feb');
save( 'C:\Users\lucak\Desktop\clash_script_upload\Second_Phase\DV\sum_nb_obj','sum_nb_obj');


% check and plot

[x y z stats]=ttest(c_feb)


% plot regression effects
model_series = mean(c_feb(:,y<0.05));
model_error = std(c_feb(:,y<0.05))./sqrt(size(c_feb,1));
para={'Value Diff','PVD','patch choice','no brainer', 'trial count','switch'};

%subplot(2,1,1)
hold on
bar(model_series','k')
set(gca,'xTick',1:numel(model_series))
set(gca,'xTickLabel',para,'FontSize',50)
set(gca,'TickLabelInterpreter','none')
%set(gca,'yTickLabel',[]);
ylabel ('Regression Weights (a.u.)')
xtickangle(45)
errorbar(1:numel(model_series),model_series,model_error,'.','color','black','linestyle','none','LineWidth',3,'CapSize',0)


ax = gca;
ax.XColor = 'k'; 
ax.YColor = 'k'; 
set(gca,'linewidth',2)
set (gca,'tickdir','out');
xlim([0.5 6.5]);
ylim([-0.7 0.3])
L = get(gca,'YLim');
set(gca,'YTick',[L(1):0.3:L(2)]);

x=0:20;

check=mes(c_feb,0,'U3_1');



