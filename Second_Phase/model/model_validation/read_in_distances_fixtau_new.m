%% This script plots Supplementary Figure S4. You will need the Gramm toolbox.
% If you have questions or comments, don`t hesitate to write an email to
% kaiserl[at]hhu.de.
close all
addpath( 'C:\Users\lucak\Desktop\fancymatplot\gramm-master')
distances_all=[]
cd('C:\Users\lucak\Desktop\clash_final\Second_Phase\Modelling_lokal\generate_recover\fixtau_multiSU_new')
vp = {'1'; '2'; '3'; '4'; '5'; '7'; '9'; '10'; '12';  '14'; '15'; '16'; '17'; '18';  '20';  '22'; '23'; '24';  '26'; '27'; '28'; '29'; '30'; '31'; '32'; '33'; '34'; '35'; '38';};



for iteration=1:numel(vp)
    data=importdata(['sim_ourdistort_fixtaufixtau_multi_' num2str(iteration) 'lu.mat']);
    
    distances_all=[distances_all;data.distances];
    
end



for num_var=1:2
    custom_map=[129 21 133;
        34 139 34]./255;
    
    figure()
    g=gramm('x',((distances_all(:,num_var))),'y',((distances_all(:,num_var+2))));
    find_border=sort(distances_all(:,num_var)-distances_all(:,num_var+2));
    border_min=find_border(0.10*length(find_border));
    border_max=find_border(0.90*length(find_border));
    diff=border_max-border_min;step=diff/5;    
    
    g.set_color_options('hue_range',[-60 60],'chroma',40,'lightness',90);    
    g.set_names('x','Real Parameters','y','Recovered Parameters');
    g.set_color_options('map',custom_map(num_var,:));
    
    
    corri=corr(distances_all(:,num_var),distances_all(:,num_var+2));
    
    g.set_title(['Recovery' num2str(num_var) '_' num2str(round(corri,2))]);
    if num_var==2 || num_var==3 || num_var==1
        g.stat_cornerhist('edges',-2.5:0.1:2.5,'aspect',0.6);
    else
        g.stat_cornerhist('edges',-5:0.1:5,'aspect',0.6);
    end
    g.geom_abline();
    g.geom_point('alpha',0.2);
    g.set_point_options('base_size',3);
    g.set_text_options('font','Arial','base_size',18);
    figure('Position',[100 100 800 600]);
    
    
    g.draw()
    
       
    if num_var==2 || num_var==3 || num_var==1
        set([g.results.stat_cornerhist.child_axe_handle],'XTick',[-2 -1 0 1 2],'fontsize',12)
    end
    
    
    
    
    
    % g.export('file_name',['gramm_export' num2str(num_var)],'file_type','pdf')
    % g.export('file_name',['gramm_export' num2str(num_var)],'file_type','jpg')
    
end

%% plot correlations between recovered params


custom_map=[0 206 209]/255


figure()
g=gramm('x',((distances_all(:,3))),'y',((distances_all(:,4))));

g.set_color_options('hue_range',[-60 60],'chroma',40,'lightness',90);
g.set_names('x','alpha','y','gamma');
g.set_color_options('map',custom_map);


corri=corr(distances_all(:,3),distances_all(:,4));

g.set_title(['Correlation recovered parameters ' '_' num2str(round(corri,2))]);
%g.stat_glm();
g.geom_abline();
g.geom_point('alpha',0.2);
g.set_point_options('base_size',3);
g.set_text_options('font','Arial','base_size',18);
figure('Position',[100 100 800 600]);


g.draw()

%  plot(g.results.stat_cornerhist.child_axe_handle,[-2 2],[0 50],'k:','LineWidth',2)
%g.export('file_name',['gramm_export_cor_variables'],'file_type','pdf')
%g.export('file_name',['gramm_export_cor_variables'],'file_type','jpg')


