%% Find out in how many trials the participant chose the correct option
% For questions or criticism please contact kaiserl@hhu.de.

%==========================================================================

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

clear all; close all; clc;

vp = {'1'; '2'; '3'; '4'; '5'; '7'; '9'; '10'; '12';  '14'; '15'; '16'; '17'; '18';  '20';  '22'; '23'; '24';  '26'; '27'; '28'; '29'; '30'; '31'; '32'; '33'; '34'; '35'; '38'};

ev_all=[];
for iteration=1:numel(vp)
    
    data=importdata(['C:\Users\lucak\Desktop\clash_script_upload\logfiles\clashMEG',vp{iteration},'.txt']);
    data=data.data;
    data=data(data(:,5)~=999,:);     % get rid of pauses
    ntrl=size(data,1);
    
    rp=data(:,[10 12])/10;
    rm=data(:,[9 11]);
    
    choices=data(:,13);
    rt=data(:,14);
    
    ev=rp.*rm;
    
    % in how many trials is one expected value better than the other one?    
    all=ev(:,1)~=ev(:,2);
    
    ev_red=ev(all,:); % only rows that have an actual better value
    maxev=double(ev_red(:,1)>ev_red(:,2));
    maxev(ev_red(:,1)<ev_red(:,2))=2;
    
    
    objcorr=sum(choices(all)==maxev); %how many times do choices match the better value?
    objcorP(iteration,:) = 100*objcorr/sum(all); %get the percentage of ev-better choices!
    
    
    clear  smaxall_grid smax_red_grid smax_red_pso smaxev_grid ev_red maxevobjcorr subcorr_pso subcorr_grid
    
end

save( 'C:\Users\lucak\Desktop\clash_script_upload\Second_Phase\DV\objcorP','objcorP');



