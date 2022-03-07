%% Sample response processing
clc
clear all;
subjects=1:14;
Blocks={'1_','2_','3_','4_','5_','6_','7_','8_','9_','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30'};
load('Order_of_blocks_for_Psych.mat');

accuracy_att=nan*ones(length(Blocks),length(subjects));
accuracy_active_att=nan*ones(length(Blocks),length(subjects));
accuracy_monit_att=nan*ones(length(Blocks),length(subjects));
accuracy_unatt=nan*ones(length(Blocks),length(subjects));
accuracy_active_unatt=nan*ones(length(Blocks),length(subjects));
accuracy_monit_unatt=nan*ones(length(Blocks),length(subjects));

TPR_att=nan*ones(length(Blocks),length(subjects));
TPR_active_att=nan*ones(length(Blocks),length(subjects));
TPR_monit_att=nan*ones(length(Blocks),length(subjects));
TPR_unatt=nan*ones(length(Blocks),length(subjects));
TPR_active_unatt=nan*ones(length(Blocks),length(subjects));
TPR_monit_unatt=nan*ones(length(Blocks),length(subjects));


FPR_att=nan*ones(length(Blocks),length(subjects));
FPR_active_att=nan*ones(length(Blocks),length(subjects));
FPR_monit_att=nan*ones(length(Blocks),length(subjects));
FPR_unatt=nan*ones(length(Blocks),length(subjects));
FPR_active_unatt=nan*ones(length(Blocks),length(subjects));
FPR_monit_unatt=nan*ones(length(Blocks),length(subjects));

TNR_att=nan*ones(length(Blocks),length(subjects));
TNR_active_att=nan*ones(length(Blocks),length(subjects));
TNR_monit_att=nan*ones(length(Blocks),length(subjects));
TNR_unatt=nan*ones(length(Blocks),length(subjects));
TNR_active_unatt=nan*ones(length(Blocks),length(subjects));
TNR_monit_unatt=nan*ones(length(Blocks),length(subjects));

FNR_att=nan*ones(length(Blocks),length(subjects));
FNR_active_att=nan*ones(length(Blocks),length(subjects));
FNR_monit_att=nan*ones(length(Blocks),length(subjects));
FNR_unatt=nan*ones(length(Blocks),length(subjects));
FNR_active_unatt=nan*ones(length(Blocks),length(subjects));
FNR_monit_unatt=nan*ones(length(Blocks),length(subjects));

dprime_pattern_active_att=nan*ones(length(Blocks),length(subjects));
dprime_pattern_monit_att=nan*ones(length(Blocks),length(subjects));
dprime_pattern_active_unatt=nan*ones(length(Blocks),length(subjects));
dprime_pattern_monit_unatt=nan*ones(length(Blocks),length(subjects));

correct_reaction_times_att=nan(length(Blocks),length(subjects));
correct_reaction_active_att=nan*ones(length(Blocks),length(subjects));
correct_reaction_monit_att=nan*ones(length(Blocks),length(subjects));

correct_reaction_times_unatt=nan(length(Blocks),length(subjects));
correct_reaction_active_unatt=nan*ones(length(Blocks),length(subjects));
correct_reaction_monit_unatt=nan*ones(length(Blocks),length(subjects));


Active_monitoring_block=[];
Cued_color_block=[];
Response_to_dots=[];
BlockNumber=[];
Distances_traj=[];
Num_double_press=0;
Num_before_appearance_association=0;



for Subj_Num=[14]

address='N:\Eye-tracking_lab\Vigilance\Newest_no_trajectory_Regine';
dirs=dir(address);
cd(address);

if Subj_Num==5
    blk_nums=25;
elseif Subj_Num==7
    blk_nums=19;
else
    blk_nums=30;    
end


    for blk=1:blk_nums
%         clearvars -except blk_nums Distances_traj Num_before_appearance_association Num_double_press sums nums BlockNumber Response_to_dots Cued_color_block Active_monitoring_block green_red_color up_down_direct_defl_det left_right_direct_app_det left_right_direct_app First_color_blocks Num_moving_dots Number_of_trials_in_blocks dprime_pattern_active dprime_pattern_monit TNR_active FNR_active FPR_active FPR_monit FPR_top FPR_bottom TPR_active TPR_monit TPR_top TPR_bottom correct_reaction_monit correct_reaction_active correct_reaction_times_top correct_reaction_times_bottom accuracy_active accuracy_monit accuracy_top accuracy_bottom Subj_Num dirs Blocks blk accuracy_top accuracy_bottom error_top error_bottom
        for i=1:size(dirs,1)
            if Subj_Num<10
                num_chars=13;
            else
                num_chars=14;
            end
            if strncmp(dirs(i).name,['Subj_',num2str(Subj_Num),'_Blk_',Blocks{blk}],num_chars)
                load(dirs(i).name);
                break;
            end
        end
        ActMon=str2double(dirs(i).name(end-8));
        mean_sampling_time=1./60;
        for dot_num=1:Num_moving_dots*Number_of_trials_in_blocks
            tr=ceil(dot_num./Num_moving_dots);
            dot_in_trial=dot_num-(tr-1).*Num_moving_dots;
            
            if ~isempty(find(key_pressed1(dot_in_trial,:,tr),1))
                
                if length(find(key_pressed1(dot_in_trial,:,tr)))>1
                    Num_double_press=Num_double_press+1;
                end
                
                key_press_sample=find(key_pressed1(dot_in_trial,:,tr), 1, 'first');
                if isnan(distance_traj1(dot_num,key_press_sample))
                    Num_before_appearance_association=Num_before_appearance_association+1;
                    distance_traj1(dot_num,key_press_sample)=3000;
                end
                dist_relative_to_boundary(dot_in_trial,tr)=distance_traj1(dot_num,key_press_sample)-hitting_border_distance;
            else
                dist_relative_to_boundary(dot_in_trial,tr)=nan;
            end
            distance_change_per_sample(dot_in_trial,tr)=(distance_traj1(dot_num,appearance_time(dot_in_trial,tr)+10)-distance_traj1(dot_num,appearance_time(dot_in_trial,tr)+20))./(11);
            
            if ~isempty(find(key_pressed2(dot_in_trial,:,tr),1))
                
                if length(find(key_pressed2(dot_in_trial,:,tr)))>1
                    Num_double_press=Num_double_press+1;
                end
                
                key_press_sample2=find(key_pressed2(dot_in_trial,:,tr), 1, 'first' );
                if isnan(distance_traj2(dot_num,key_press_sample2))
                    Num_before_appearance_association=Num_before_appearance_association+1;
                    distance_traj2(dot_num,key_press_sample)=3000;
                end
                dist_relative_to_boundary2(dot_in_trial,tr)=distance_traj2(dot_num,key_press_sample2)-hitting_border_distance;
            else
                dist_relative_to_boundary2(dot_in_trial,tr)=nan;
            end
            distance_change_per_sample2(dot_in_trial,tr)=(distance_traj2(dot_num,appearance_time2(dot_in_trial,tr)+10)-distance_traj2(dot_num,appearance_time2(dot_in_trial,tr)+20))./(11);
        end
        
        
        distance_change_per_sample(distance_change_per_sample<0)=mean(distance_change_per_sample(distance_change_per_sample>0));
        distance_change_per_sample2(distance_change_per_sample2<0)=mean(distance_change_per_sample2(distance_change_per_sample2>0));
        
        reaction_times=((-dist_relative_to_boundary)./distance_change_per_sample).*mean_sampling_time;
        reaction_times2=((-dist_relative_to_boundary2)./distance_change_per_sample2).*mean_sampling_time;
        %% Behavioural Performance

        % attended
        tp_att=0;
        tn_att=0;
        fp_F_att=0;
        fp_S_att=0;
        fp_T_att=0;
        fn_att=0;
            
        g=0;        
        for dot_num=1:Num_moving_dots*Number_of_trials_in_blocks
            tr=ceil(dot_num./Num_moving_dots);
            dot_in_trial=dot_num-(tr-1).*Num_moving_dots;
            
            if sum(dot_in_trial==top_events(:,tr))==1 && dot_color(dot_in_trial,tr)==First_color_blocks(Subj_Num,blk)
                g=g+1;
                if isnan(reaction_times(dot_in_trial,tr)) && (top_events(tr)~=top_targets(tr))
                    tn_att=tn_att+1;    % number of non-target events with no resp;
                elseif ~isnan(reaction_times(dot_in_trial,tr)) && top_events(tr)~=top_targets(tr) && (reaction_times(dot_in_trial,tr)<0)
                    fp_F_att=fp_F_att+1;    % number of non-target events with fast resp;
                elseif ~isnan(reaction_times(dot_in_trial,tr)) && top_events(tr)~=top_targets(tr) && (reaction_times(dot_in_trial,tr)>=0)
                    fp_S_att=fp_S_att+1;    % number of non-target events with Slow resp;
                elseif ~isnan(reaction_times(dot_in_trial,tr)) && top_events(tr)==top_targets(tr) && (reaction_times(dot_in_trial,tr)<0)
                    fp_T_att=fp_T_att+1;    % number of target events with Too early resp;
                elseif isnan(reaction_times(dot_in_trial,tr)) && top_events(tr)==top_targets(tr)
                    fn_att=fn_att+1;    % number of target events with no resp;
                elseif ~isnan(reaction_times(dot_in_trial,tr)) && top_events(tr)==top_targets(tr) && reaction_times(dot_in_trial,tr)>0
                    tp_att=tp_att+1;    % number of target events with resp;
                    correct_reaction_times_att(blk,Subj_Num)=correct_reaction_times_att(blk,Subj_Num)+reaction_times(dot_in_trial,tr);
                end
            end
            
            if sum(dot_in_trial==top_events2(:,tr))==1 && dot_color2(dot_in_trial,tr)==First_color_blocks(Subj_Num,blk)
                g=g+1;
                if isnan(reaction_times2(dot_in_trial,tr)) && (top_events2(tr)~=top_targets2(tr))
                    tn_att=tn_att+1;
                elseif ~isnan(reaction_times2(dot_in_trial,tr)) && top_events2(tr)~=top_targets2(tr) && (reaction_times2(dot_in_trial,tr)<0)
                    fp_F_att=fp_F_att+1;
                elseif ~isnan(reaction_times2(dot_in_trial,tr)) && top_events2(tr)~=top_targets2(tr) && (reaction_times2(dot_in_trial,tr)>=0)
                    fp_S_att=fp_S_att+1;
                elseif ~isnan(reaction_times2(dot_in_trial,tr)) && top_events2(tr)==top_targets2(tr) && (reaction_times2(dot_in_trial,tr)<0)% || reaction_times2(dot_in_trial,tr)>time_to_touch_the_obstacle2(dot_in_trial,tr))
                    fp_T_att=fp_T_att+1;
                elseif isnan(reaction_times2(dot_in_trial,tr)) && top_events2(tr)==top_targets2(tr)
                    fn_att=fn_att+1;
                elseif ~isnan(reaction_times2(dot_in_trial,tr)) && top_events2(tr)==top_targets2(tr) && reaction_times2(dot_in_trial,tr)>0 %&& reaction_times2(dot_in_trial,tr)<time_to_touch_the_obstacle2(dot_in_trial,tr)
                    tp_att=tp_att+1;
                    correct_reaction_times_att(blk,Subj_Num)=correct_reaction_times_att(blk,Subj_Num)+reaction_times2(dot_in_trial,tr);
                end
            end
        end
        
        fp_att=fp_F_att+fp_S_att+fp_T_att;
        sum([tp_att tn_att fp_att fn_att])

        accuracy_att(blk,Subj_Num)=(tp_att+tn_att)./(sum(top_events>0)+sum(top_events2>0));
        TPR_att(blk,Subj_Num)=(tp_att)./(tp_att+fn_att);
        TNR_att(blk,Subj_Num)=(tn_att)./(tn_att+fp_att);
        FPR_att(blk,Subj_Num)=(fp_att)./(fp_att+tn_att);
        FNR_att(blk,Subj_Num)=(fn_att)./(tp_att+fn_att);
        correct_reaction_times_att(blk,Subj_Num)=correct_reaction_times_att(blk,Subj_Num)./tp_att;
        if ActMon==1
            accuracy_active_att(blk,Subj_Num)=accuracy_att(blk,Subj_Num);
            TPR_active_att(blk,Subj_Num)=TPR_att(blk,Subj_Num);
            TNR_active_att(blk,Subj_Num)=TNR_att(blk,Subj_Num);
            FPR_active_att(blk,Subj_Num)=FPR_att(blk,Subj_Num);
            FNR_active_att(blk,Subj_Num)=FNR_att(blk,Subj_Num);
            dprime_pattern_active_att(blk,Subj_Num)=TPR_active_att(blk,Subj_Num)-FPR_active_att(blk,Subj_Num);
            correct_reaction_active_att(blk,Subj_Num)=correct_reaction_times_att(blk,Subj_Num);
        else
            accuracy_monit_att(blk,Subj_Num)=accuracy_att(blk,Subj_Num);
            TPR_monit_att(blk,Subj_Num)=TPR_att(blk,Subj_Num);
            TNR_monit_att(blk,Subj_Num)=TNR_att(blk,Subj_Num);
            FPR_monit_att(blk,Subj_Num)=FPR_att(blk,Subj_Num);
            FNR_monit_att(blk,Subj_Num)=FNR_att(blk,Subj_Num);
            dprime_pattern_monit_att(blk,Subj_Num)=TPR_monit_att(blk,Subj_Num)-FPR_monit_att(blk,Subj_Num);
            correct_reaction_monit_att(blk,Subj_Num)=correct_reaction_times_att(blk,Subj_Num);
        end
        
% unattended
        tp_unatt=0;
        fp_unatt=0;
        
        for dot_num=1:Num_moving_dots*Number_of_trials_in_blocks
            tr=ceil(dot_num./Num_moving_dots);
            dot_in_trial=dot_num-(tr-1).*Num_moving_dots;
            if sum(dot_in_trial==top_events(:,tr))==1 && dot_color(dot_in_trial,tr)~=First_color_blocks(Subj_Num,blk)
                if isnan(reaction_times(dot_in_trial,tr))
                    tp_unatt=tp_unatt+1;    % number of non-target events with no resp;
                else
                    fp_unatt=fp_unatt+1;    % number of target events with resp;
                    correct_reaction_times_unatt(blk,Subj_Num)=correct_reaction_times_unatt(blk,Subj_Num)+reaction_times(dot_in_trial,tr);
                end
            end
            
            if sum(dot_in_trial==top_events2(:,tr))==1 && dot_color2(dot_in_trial,tr)~=First_color_blocks(Subj_Num,blk)
                if isnan(reaction_times2(dot_in_trial,tr))
                    tp_unatt=tp_unatt+1;
                else
                    fp_unatt=fp_unatt+1;
                    correct_reaction_times_unatt(blk,Subj_Num)=correct_reaction_times_unatt(blk,Subj_Num)+reaction_times2(dot_in_trial,tr);
                end
            end
        end
        
        
%         accuracy_unatt(blk,Subj_Num)=(tp_unatt+tn_unatt)./(sum(top_events>0)+sum(top_events2>0));
        TPR_unatt(blk,Subj_Num)=(tp_unatt)./(tp_unatt+fp_unatt);
%         TNR_unatt(blk,Subj_Num)=(tn_unatt)./(tn_unatt+fp_unatt);
        FPR_unatt(blk,Subj_Num)=(fp_unatt)./(fp_unatt+tp_unatt);
%         FNR_unatt(blk,Subj_Num)=(fn_unatt)./(tp_unatt+fn_unatt);
        correct_reaction_times_unatt(blk,Subj_Num)=correct_reaction_times_unatt(blk,Subj_Num)./tp_unatt;
        
        if ActMon==1
%             accuracy_active_unatt(blk,Subj_Num)=accuracy_unatt(blk,Subj_Num);
            TPR_active_unatt(blk,Subj_Num)=TPR_unatt(blk,Subj_Num);
            TNR_active_unatt(blk,Subj_Num)=TNR_unatt(blk,Subj_Num);
            FPR_active_unatt(blk,Subj_Num)=FPR_unatt(blk,Subj_Num);
            FNR_active_unatt(blk,Subj_Num)=FNR_unatt(blk,Subj_Num);
            dprime_pattern_active_unatt(blk,Subj_Num)=TPR_active_unatt(blk,Subj_Num)-FPR_active_unatt(blk,Subj_Num);
            correct_reaction_active_unatt(blk,Subj_Num)=correct_reaction_times_unatt(blk,Subj_Num);
        else
%             accuracy_monit_unatt(blk,Subj_Num)=accuracy_unatt(blk,Subj_Num);
            TPR_monit_unatt(blk,Subj_Num)=TPR_unatt(blk,Subj_Num);
            TNR_monit_unatt(blk,Subj_Num)=TNR_unatt(blk,Subj_Num);
            FPR_monit_unatt(blk,Subj_Num)=FPR_unatt(blk,Subj_Num);
            FNR_monit_unatt(blk,Subj_Num)=FNR_unatt(blk,Subj_Num);
            dprime_pattern_monit_unatt(blk,Subj_Num)=TPR_monit_unatt(blk,Subj_Num)-FPR_monit_unatt(blk,Subj_Num);
            correct_reaction_monit_unatt(blk,Subj_Num)=correct_reaction_times_unatt(blk,Subj_Num);
        end
        
    end
    
%     left_right_direct_app=reshape(left_right_direct_app',[size(left_right_direct_app,1)*size(left_right_direct_app,2) 1]);
%     left_right_direct_app_det=reshape(left_right_direct_app_det',[size(left_right_direct_app_det,1)*size(left_right_direct_app_det,2) 1]);
%     up_down_direct_defl_det=reshape(up_down_direct_defl_det',[size(up_down_direct_defl_det,1)*size(up_down_direct_defl_det,2) 1]);
%     green_red_color=reshape(green_red_color',[size(green_red_color,1)*size(green_red_color,2) 1]);
%     save(['Event_labels_Subj_',num2str(Subj_Num),'.mat'],'Distances_traj','BlockNumber','Response_to_dots','Active_monitoring_block','Cued_color_block','left_right_direct_app','left_right_direct_app_det','up_down_direct_defl_det','green_red_color')
end
[Num_double_press Num_before_appearance_association]
save('Beh_data_01_Subjs_separated_att_unat.mat','TPR_active_att','TNR_active_att','FPR_active_att','FNR_active_att','TPR_monit_att','TNR_monit_att','FPR_monit_att','FNR_monit_att','correct_reaction_active_att','correct_reaction_monit_att','TPR_active_unatt','TNR_active_unatt','FPR_active_unatt','FNR_active_unatt','TPR_monit_unatt','FPR_monit_unatt','correct_reaction_active_unatt','correct_reaction_monit_unatt')
ccc
%% Behav data analysis: pooling first and second halves for the subject
% load('Beh_data_14_Subjs.mat');
% figure;
% dprime_active=dprime(nanmean(TPR_active,2),nanmean(FPR_active,2));
% plot(mean([dprime_active(1:15) dprime_active(16:end)],2),'LineWidth',1.5)
% hold on;
% dprime_monit=dprime(nanmean(TPR_monit,2),nanmean(FPR_monit,2));
% plot(mean([dprime_monit(1:15) dprime_monit(16:end)],2),'LineWidth',1.5)
% grid on;
% xlabel ('Block number');
% ylabel ('dprime')
% 
% figure;
% TPR_active=nanmean(TPR_active,2);
% plot(mean([TPR_active(1:15) TPR_active(16:end)],2),'LineWidth',1.5)
% shadedErrorBar([1:15],TPR_active(1:15),{@median,@std},{'r-o','markerfacecolor','r'})
% hold on;
% TPR_monit=nanmean(TPR_monit,2);
% plot(mean([TPR_monit(1:15) TPR_monit(16:end)],2),'LineWidth',1.5)
% % plot(nanmean([mean([TPR_active(1:15) TPR_active(16:end)],2) mean([TPR_monit(1:15) TPR_monit(16:end)],2)],2),'LineWidth',1.5)
% grid on;
% % legend Active Monit
% xlabel ('Block number');
% ylabel ('TPR')
% 
% figure;
% FPR_active=nanmean(FPR_active,2);
% plot(mean([FPR_active(1:15) FPR_active(16:end)],2),'LineWidth',1.5)
% hold on;
% FPR_monit=nanmean(FPR_monit,2);
% plot(mean([FPR_monit(1:15) FPR_monit(16:end)],2),'LineWidth',1.5)
% % plot(nanmean([mean([FPR_active(1:15) FPR_active(16:end)],2) mean([FPR_monit(1:15) FPR_monit(16:end)],2)],2),'LineWidth',1.5)
% grid on;
% % legend Active Monit
% xlabel ('Block number');
% ylabel ('FPR')
% 
% figure;
% correct_reaction_active=nanmean(correct_reaction_active,2);
% plot(nanmean([correct_reaction_active(1:15) correct_reaction_active(16:end)],2),'LineWidth',1.5)
% hold on;
% correct_reaction_monit=nanmean(correct_reaction_monit,2);
% plot(nanmean([correct_reaction_monit(1:15) correct_reaction_monit(16:end)],2),'LineWidth',1.5)
% % plot(nanmean([mean([correct_reaction_active(1:15) correct_reaction_active(16:end)],2) mean([correct_reaction_monit(1:15) correct_reaction_monit(16:end)],2)],2),'LineWidth',1.5)
% grid on;
% % legend Active Monit
% xlabel ('Block number');
% ylabel ('Reaction time [s]')

%%
clc;
clear all;
% close all;
load('Beh_data_01_Subjs_separated_att_unat.mat');
RT=0;
if RT==0
    dataA1=TPR_active_att;
    dataB1=TPR_monit_att;
else   
    dataA1=correct_reaction_active_att;
    dataB1=correct_reaction_monit_att;
end

MeanA=(nanmean([dataA1(1:15,:) dataA1(16:30,:)],2));
StdA=nanstd([dataA1(1:15,:) dataA1(16:30,:)]');
MeanM=(nanmean([dataB1(1:15,:) dataB1(16:30,:)],2));
StdM=nanstd([dataB1(1:15,:) dataB1(16:30,:)]');
figure;
Shad1=shadedErrorBar([1:15],MeanA,(StdA)./sqrt(size(dataA1,2))*1.96,{'color',[0.1 0.1 0.9],'LineWidth',3},1);
hold on;
Shad2=shadedErrorBar([1:15],MeanM,(StdA)./sqrt(size(dataA1,2))*1.96,{'color',[0.9 0.1 0.1],'LineWidth',3},1);
for block=1:15
    tmpA=[dataA1(1:15,:) dataA1(16:30,:)];
    tmpB=[dataB1(1:15,:) dataB1(16:30,:)];
    Effects(block)= bf.ttest2(tmpA(block,~isnan(tmpA(block,:)))',tmpB(block,~isnan(tmpB(block,:)))');
    if Effects(block)>10
        Bayes(block)=2.5;
    elseif Effects(block)>3 && Effects(block)<=10
        Bayes(block)=1.5;
    elseif Effects(block)>1 && Effects(block)<=3
        Bayes(block)=0.5;
    elseif Effects(block)<1 && Effects(block)>=1/3
        Bayes(block)=-0.5;
    elseif Effects(block)<1/3 && Effects(block)>=1/10
        Bayes(block)=-1.5;
    elseif Effects(block)<1/10
        Bayes(block)=-2.5;
    end
end

if RT
    Baseline=0.29;
    steps=0.0015;
    distans=1.5; % times step
    ylabel('Reaction time (s)')
    ylim([0.283 0.4])

else
    Baseline=0;
    steps=0.01;
    distans=1.5; % times step
    ylabel('Proportion of misses')
    ylim([-0.05 0.8])
end

for block=1:size(Bayes,2)
    if Bayes(block)==-0.5 || Bayes(block)==0.5
        plots(block)=plot(block,Bayes(block).*steps+Baseline,'LineStyle','none','marker','o','Color','k','linewidth',2,'markersize',10);
    else
        plots(block)=plot(block,Bayes(block).*steps+Baseline,'LineStyle','none','marker','o','MarkerFaceColor','k','Color','k','linewidth',2,'markersize',10);
    end
    hold on;
end
baseline_temp=Baseline;
line([1 length(Bayes)],[baseline_temp baseline_temp],'linestyle','-.','Color','k','linewidth',2);
line([1 length(Bayes)],[baseline_temp baseline_temp]-steps,'Color','k','linewidth',2);
line([1 length(Bayes)],[baseline_temp baseline_temp]-2*steps,'Color','k','linewidth',2);
line([1 length(Bayes)],[baseline_temp baseline_temp]-3*steps,'Color','k','linewidth',2);
line([1 length(Bayes)],[baseline_temp baseline_temp]+steps,'Color','k','linewidth',2);
line([1 length(Bayes)],[baseline_temp baseline_temp]+2*steps,'Color','k','linewidth',2);
line([1 length(Bayes)],[baseline_temp baseline_temp]+3*steps,'Color','k','linewidth',2);

xlabel('Block #')
set(gca,'FontSize',30);
xlim([0 16])

