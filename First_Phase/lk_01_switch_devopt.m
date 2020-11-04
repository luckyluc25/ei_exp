%% Are there differences in patch value difference as a function of different cost levels?

% For questions, comments or criticism, please contact kaiserl@hhu.de

% required input: behavioral data (.txt-file with column headers), columns:
% 1 =   subject ID
% 2 =   trial no
% 3 =   value of patch A
% 4 =   value of patch B
% 5 =   position of patches on screen (1=A left, B right, 2=A right, B left) Caution: 999=Pause
% 6 =   patch choice on current trial
% 7 =   switch cost
% 8 =   RT for patch decision
% 9 =  2nd stage: reward magnitude left
% 10=  2nd stage: reward probability left
% 11=  2nd stage: reward magnitude right
% 12=  2nd stage: reward probability right
% 13=  choice (left/right) at 2nd stage
% 14=  RT for 2nd stage choice
% 15=  outcome (reward/no reward) left for 2nd stage
% 16=  outcome (reward/no reward) right for 2nd stage
% 17=  cumulative points won by subject


%last edited 9/7/20

%% TO DO BEFORE RUNNING THE SCRIPT
%  Depends on the MEST toolbox --> download at: https://de.mathworks.com/matlabcentral/fileexchange/32398-hhentschke-measures-of-effect-size-toolbox
%  change your paths accordingly

%% %% ==================================================================%


function[]=lk_01_switch_devopt()

%clear up before you start
clear all; close all;
addpath('C:\Users\luca\Desktop\magma\MatFunc');
addpath('C:\Users\lucak\Desktop\measures-of-effect-size-toolbox-master');

%Load data of participants with valid MRS measurements in all voxels of interest-be careful MRS 6 is MEG and behav 7!
vp = {'1'; '2'; '3'; '4'; '5'; '7'; '9'; '10'; '12';  '14'; '15'; '16'; '17'; '18';  '20';  '22'; '23'; '24';  '26'; '27'; '28'; '29'; '30'; '31'; '32'; '33'; '34'; '35'; '38'};

%loop through all subjects
for iteration=1:numel(vp)
    
    data=importdata(['C:\Users\lucak\Desktop\clash_script_upload\logfiles\clashMEG',vp{iteration},'.txt']);
    data=data.data;
    data=data(data(:,5)~=999,:);     % get rid of pauses
    ntrl=size(data,1);
    
    sw = [0; abs(diff(data(:,6)))];  % outcome variable: switch or stay in patch on current trial
    
    sum_switches (iteration,:)=sum(sw); switches_all_per(iteration,:)=sum(sw)/size(sw,1);% how many times did participants switch?
    
    vcurr=zeros(ntrl,1);             %We use patch value difference from previous trials since subjects have not seen updated pvd at time of choice
    valt=zeros(ntrl,1);
    
    for k=2:ntrl
        if data(k-1,6)==1
            vcurr(k) = data(k-1,3);
            valt(k) = data(k-1,4);
        elseif data(k-1,6)==2
            vcurr(k) = data(k-1,4);
            valt(k) = data(k-1,3);
        end
    end
    
    
    pvd = valt-vcurr;
    c=data(:,7);                     %c=switch costs
    
    
    % find switches for each cost level
    
    if ~isempty(find(sw==1 & c==5))
        vd1=pvd(find(sw==1 & c==5));
    else
        vd1=NaN;
    end
    if ~isempty(find(sw==1 & c==10))
        vd2=pvd(find(sw==1 & c==10));
    else
        vd2=NaN;
    end
    if ~isempty(find(sw==1 & c==15))
        vd3=pvd(find(sw==1 & c==15));
    else
        vd3=NaN;
    end
    if ~isempty(find(sw==1 & c==20))
        vd4=pvd(find(sw==1 & c==20));
    else
        vd4=NaN;
    end
    
    delta_med(iteration,:)= [nanmedian(repmat(vd1,2,1))  nanmedian(repmat(vd2,2,1))    nanmedian(repmat(vd3,2,1))    nanmedian(repmat(vd4,2,1))];
    
    
end

%% save and show results

save(  'C:\Users\lucak\Desktop\clash_script_upload\First_Phase\DV\sum_switches.mat', 'sum_switches');
save(  'C:\Users\lucak\Desktop\clash_script_upload\First_Phase\DV\switches_all_per', 'switches_all_per');
save(  'C:\Users\lucak\Desktop\clash_script_upload\First_Phase\DV\delta_med', 'delta_med');


%% How often did subjects switch?
avg_switches=mean(switches_all_per); std_avg_switches=std(switches_all_per)./sqrt(size(switches_all_per,1));
disp(['Mean per switches is', num2str(avg_switches), '_',num2str(std_avg_switches)]);

%% Is the median patch value difference cost dependent?
clear delta
[row,col]=find(isnan(delta_med(:,:)));
delta_med(row,:)=[];

mean_val_diff=mean(delta_med);sem_val_diff=std(delta_med)./sqrt(size(delta_med,1));
disp(['Mean val diff is', num2str(mean_val_diff), '_',num2str(sem_val_diff)]);

%plot results
figure;
hbar=bar(nanmean(delta_med),'k'); hold all;
set(gca,'tickdir','out','box','off','fontsize',65,'XTickLabel',{'5','10','15','20'});
xlabel('cost'),ylabel ('average value difference');
ylim([27 42]);
hold on
errorbar(1:4,nanmean(delta_med),nanstd(delta_med)./sqrt(size(delta_med,1)),'k','linestyle','none','linewidth',3,'CapSize',0)
ax = gca;
ax.XColor = 'k'; 
ax.YColor = 'k'; 
set(gca,'linewidth',2)
xlim([0.5 4.5])

disp (mean(delta_med)); disp (std(delta_med)./sqrt(size(delta_med,1)))

% run RM-ANOVA
t = table((1:length(delta_med))',delta_med(:,1),delta_med(:,2),delta_med(:,3),delta_med(:,4),'VariableNames',{'Subjects','cost1','cost2','cost3','cost4'});
Costs = table([1 2 3 4]','VariableNames',{'Costs'});
rm = fitrm(t,'cost1-cost4~1','WithinDesign',Costs);
ranovatbl = ranova(rm);

data = [  delta_med(:)        squash(repmat(1:size(delta_med,2),size(delta_med,1),1))  repmat((1:size(delta_med,1))',size(delta_med,2),1)];
RMANOVA1(data)


% fit linear trend
des = [ones(4,1) (1:4)'];
copes=pinv(des)*delta_med'; copes=copes';

[h,p,ci,stats] = ttest(copes);

%test effect sizes
check=mes(copes,0,'U3_1');
delta_med_all=delta_med(:);group=[repmat(1,length(delta_med),1);repmat(2,length(delta_med),1);repmat(3,length(delta_med),1);repmat(4,length(delta_med),1)];
mes1way(delta_med_all,'eta2','group',group,'isDep',1)

end






