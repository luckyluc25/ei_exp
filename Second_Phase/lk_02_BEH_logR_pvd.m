
%% Do a logistic regression on choice to estimate which factors influence whether participants choose the right option
%for questions please contact kaiserl@hhu.de

%=========================================================================
% required input: behavioral data (.txt-file with column headers), columns have are organized:
% 1 =   subject ID
% 2 =   trial no
% 3 =   value of patch A
% 4 =   value of patch B
% 5 =   position of patches on screen (1=A left, B right, 2=A right, B left) Caution: 999=Pause
% 6 =   patch choice on current trial
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

% add the path you need
addpath('C:\Users\lucak\Desktop\MatFunc');
addpath('C:\Users\luca\Desktop\hhentschke-measures-of-effect-size-toolbox-3d90ae5')
addpath( 'C:\Users\lucak\Desktop\measures-of-effect-size-toolbox-master')

clear all; close all; clc;

vp = {'1'; '2'; '3'; '4'; '5'; '7'; '9'; '10'; '12';  '14'; '15'; '16'; '17'; '18';  '20';  '22'; '23'; '24';  '26'; '27'; '28'; '29'; '30'; '31'; '32'; '33'; '34'; '35'; '38'};

for iteration=1:numel(vp)
    
    data=importdata(['C:\Users\lucak\Desktop\clash_script_upload\logfiles\clashMEG',vp{iteration},'.txt']);
    data=data.data;
    data=data(data(:,5)~=999,:);     % get rid of pauses
    ntrl=size(data,1);
    
    
    rp=data(:,[10 12])/10;
    rm=data(:,[9 11]);
    choices=data(:,13); choices(choices==1) = 0;
    choices(choices==2) = 1;
    
    sw = [0; abs(diff(data(:,6)))];
    
    ev=rp.*rm;
    
    % extract pvd choices
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
    
    test=unique(choices_pvd);
    
    if numel(test)~=2
        keyboard
    end
    
    
    rt=data(:,14);
    
    
    index=3;
    choice_right=double(choices==1);
    repchoices= toeplitz(choice_right',zeros(1,index+1));
    repchoices=repchoices(:,2:end);
    
    
    all_wins=double((data(:,16)==1 & data(:,13)==2) | (data(:,15)==1 & data(:,13)==1));
    all_wins=toeplitz(all_wins,zeros(1,index+1));
    all_wins=all_wins(:,2:end);
    
    % we take the current patch value here since they are already known by
    % the subject during patch leaving decisions!
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
    nb_obj = double((rp(:,1)>=rp(:,2) & rm(:,1)>=rm(:,2)) | (rp(:,2)>=rp(:,1) & rm(:,2)>=rm(:,1)));
    trial_count=[1:length(nb_obj)]';
    
    
    
    x=zscore([(ev(:,2)-ev(:,1)) ev(:,1)+ev(:,2) costs (pvd') repchoices(:,1) all_wins(:,1) choices_pvd nb_obj trial_count sw]);
    
    
    y = (choices==1);
    c = glmfit(x,y,'binomial','link','logit');
    
    
    logRes_pvd(iteration,:) = (c);
    
    clear x y rm rp
end

% save and plot
save( 'C:\Users\lucak\Desktop\clash_script_upload\Second_Phase\DV\logRes_pvd','logRes_pvd');


mean(logRes_pvd); std(logRes_pvd./sqrt(size(logRes_pvd,1)))
[x y z stats]=ttest(logRes_pvd)
check=mes(logRes_pvd,0,'U3_1')

%% plot effects on choice

model_series = mean(logRes_pvd);
model_error = std(logRes_pvd)./sqrt(size(logRes_pvd,1));

model_series_short = mean(logRes_pvd(:,y<0.05));
model_error_short = std(logRes_pvd(:,y<0.05))./sqrt(size(logRes_pvd(:,y<0.05),1));

hold on
bar(model_series_short','k')
set(gca,'xTick',1:numel(model_series_short))
set(gca,'xTickLabel',para,'FontSize',50)
set(gca,'TickLabelInterpreter','none')
ylabel ('Regression Weights (a.u.)')
errorbar(1:numel(model_series_short),model_series_short,model_error_short,'.','color','black','linestyle','none','LineWidth',3,'CapSize',0)
ax = gca;
ax.XColor = 'k'; 
ax.YColor = 'k';
set(gca,'linewidth',2,'tickdir','out')
ylim([-0.2 1.2]);
xlim([0.5 4.5]);
L=get(gca,'YLim');
set(gca,'YTick',[L(1):0.2:L(2)])


