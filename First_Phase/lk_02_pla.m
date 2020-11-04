%% find out subject specific differences in 'patch leaving advantages'

% required input: behavioral data (.txt-file with column headers), columns:
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

%==================================================================%

%% Please contact kaiserl@hhu.de for questions or comments

%last edited 5/3/20

%Please change paths accordingly before running the script.

% If you have questions or criticism, please contact kaiserl@hhu.de
%==================================================================%

clear all; close all;
addpath( 'C:\Users\luca\Desktop\MatFunc');

vp = {'1'; '2'; '3'; '4'; '5'; '7'; '9'; '10'; '12';  '14'; '15'; '16'; '17'; '18';  '20';  '22'; '23'; '24';  '26'; '27'; '28'; '29'; '30'; '31'; '32'; '33'; '34'; '35'; '38'};


for iteration=1:numel(vp)
    
    data=importdata(['C:\Users\lucak\Desktop\clash_script_upload\logfiles\clashMEG',vp{iteration},'.txt']);
    
    data=data.data;
    data=data(data(:,5)~=999,:);     % get rid of pauses
    ntrl=size(data,1);
    
    %find out how many times subjects switched patches
    sw = [0; abs(diff(data(:,6)))];          % outcome variable: switch or stay in patch on current trial
    sw_plot(iteration,:)=sw;
    numel_switches(iteration)=sum(sw);
    
    %find out patch value differences that made participants leave their patch
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
    
    costs=data(:,7); pvd=valt-vcurr;
    
    
    %patch leaving indices depend on patch value differences and costs
    pla=[pvd costs];
    distance=pla(sw==1,:);
    
    switch_index(iteration,:)=mean(distance(:,1)-distance(:,2));
    
    distance_to_plot{iteration,:}=pla(:,1)-pla(:,2);
    
    
end

save( 'C:\Users\lucak\Desktop\clash_script_upload\First_Phase\DV\switch_index','switch_index');

%% plot timecourse of PLA

scatter(find(sw_plot(1,:)==1),zeros(sum(sw_plot(1,:)),1)+47,200,[0 128/256 0],'filled');
hold on
plot(distance_to_plot{1,1},'LineWidth',3,'Color','black');
set(gca,'FontSize',50)
hold on
plot  ([1 320], [0 0],'Color','black','LineWidth',3)
xlabel('trial number')
ylabel('T_o_p_t')
ax = gca;
ax.XColor = 'k'; % Red
ax.YColor = 'k'; % Blue
set(gca,'linewidth',2)
set(gca,'xlim',[0 320],'tickdir','out');
set(gca,'ylim',[-75 55]);
L = get(gca,'XLim');
set(gca,'XTick',[L(1):50:L(2)])



