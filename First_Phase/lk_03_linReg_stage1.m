%% find out subject specific patch leaving indices

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

%==================================================================%

%Please contact kaiserl@hhu.de for questions

%last edited 9/3/19

%Please change paths accordingly before running the script.

%===================================================================%

addpath('C:\Users\luca\Desktop\magma\MatFunc');
addpath ('C:\Users\lucak\Desktop\measures-of-effect-size-toolbox-master')

clear all; close all; clc;

vp = {'1'; '2'; '3'; '4'; '5'; '7'; '9'; '10'; '12';  '14'; '15'; '16'; '17'; '18';  '20';  '22'; '23'; '24';  '26'; '27'; '28'; '29'; '30'; '31'; '32'; '33'; '34'; '35'; '38'};



for iteration=1:numel(vp)
    data=importdata(['C:\Users\lucak\Desktop\clash_script_upload\logfiles\clashMEG',vp{iteration},'.txt']);
    data=data.data;
    data=data(data(:,5)~=999,:);     % get rid of pauses
    sw = [0; abs(diff(data(:,6)))];
    ntrl=size(data,1);
    
    ch_side=data(:,5); ch_side=[0; abs(diff(ch_side))];
    
    vcurr=zeros(ntrl,1);
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
    
    rt=data(:,8);
    rt_1st(:,iteration)=median(rt);
    costs=data(:,7);
    
    
    trial_count=1:length(sw);
    
    %% include previous 2nd stage parameter
    
    wins=zeros(size(rt));
    wins( (data(:,13)==1 & data(:,15)==1) | (data(:,13)==2 & data(:,16)==1) ) = 1;       % did subject win at 2nd stage?
    
    k=3;
    warning('off');
    repwins=toeplitz(wins,zeros(1,k+1));
    warning('on');
    
    repwins=repwins(:,2:end);
    
    %% do linear regression on log RT
    
    % we go on from the second trial since we always use the alternative value that has been shown in the trial before (patch values of the
    % current trial are unknown at the time of choice.
    
    
    xrt=[zscore([valt(2:end)-vcurr(2:end)  sw(2:end) (costs(2:end)) trial_count(2:end)'  ch_side(2:end) repwins(2:end,1)])];
    
    corr_check=corr(xrt);
    
    check_coeffs=fitlm(xrt,zscore(log(rt(2:end))));
    
    
    linReg_1stphase_check(iteration,:)=(check_coeffs.Coefficients.Estimate); 
    
    
    
end

save( 'C:\Users\lucak\Desktop\clash_script_upload\First_Phase\DV\linReg_1stphase_check','linReg_1stphase_check')
save( 'C:\Users\lucak\Desktop\clash_script_upload\First_Phase\DV\rt_1st','rt_1st')

%% test coefficients and plot!

[x y z stats]=ttest(linReg_1stphase_check);
%plot
model_series = mean(linReg_1stphase_check(:,2:end));
model_error = std(linReg_1stphase_check(:,2:end))./sqrt(size(linReg_1stphase_check(:,2:end),1));
para={'PVD','SW','Costs', 'trial count', 'CH','repwins'};

hold on
bar(model_series','k')
set(gca,'xTick',1:numel(model_series))
set(gca,'xTickLabel',para,'FontSize',50)
set(gca,'TickLabelInterpreter','none')
ylabel ('Regression Weights (a.u.)')
xtickangle(45)
errorbar(1:numel(model_series),model_series,model_error,'.','color','black','linestyle','none','LineWidth',3,'CapSize',0)
ax = gca;
ax.XColor = 'k';
ax.YColor = 'k';
set(gca,'linewidth',2)
set(gca,'tickdir','out')
xlim([0.5 6.5]);
ylim([-0.2 0.6])


L = get(gca,'YLim');
set(gca,'YTick',[L(1):0.2:L(2)])


%% check for effect sizes with the MES toolbox
check=mes(linReg_1stphase_check,0,'U3_1');


