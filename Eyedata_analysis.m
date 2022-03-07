clc;
% clear all;
close all;
load('Order_of_blocks_for_Psych.mat');

% Subject=3;
for Subject=[1:10 11:14 16:19]
% for Subject=[12]
    Tr=0;
    if Tr==1
        traintest='Tr';
    else
        traintest='Ts';
    end
    
    eyedata=Edf2Mat(['Subj',num2str(Subject),traintest,'.edf']);
    Xtemp=eyedata.Samples.posX;
    Ytemp=eyedata.Samples.posY;

    %% Extracting and plotting the data manually
    % Finding block times
    Data_collection=nan*ones(2,26*120000,2);
    LabelsAM=nan*ones(1,26*120000);
    LabelsTB=nan*ones(1,26*120000);
    for Block_to_analize=1:26
        clearvars -except eyedata Xtemp Ytemp Subject Block_to_analize Data_collection Pupil_sizes Top_blocks Active_blocks LabelsAM LabelsTB
        if Active_blocks(Subject,Block_to_analize)==1
            Active_Monitoring='Active';
        else
            Active_Monitoring='Monitoring';
        end
        if Top_blocks(Subject,Block_to_analize)==1
            Cued_Location='Top';
        else
            Cued_Location='Bottom';
        end
        
        block=0;
        for count=1:length(eyedata.Events.Messages.info)
            if strcmp(eyedata.Events.Messages.info{1,count},'Block_Onset')
                block=block+1;
                BlockTimes(block,1)=eyedata.Events.Messages.time(count);
            end
            if strcmp(eyedata.Events.Messages.info{1,count},'Start_of_Rest_Time')
                BlockTimes(block,2)=eyedata.Events.Messages.time(count);
            end
        end
        
        % Removing artifactual data
        % Eblink
        % Esacc
        % Efix
        desired_block_fixations=zeros(length(Xtemp),1);
        desired_block_fixations(find(eyedata.Samples.time==BlockTimes(Block_to_analize,1)):find(eyedata.Samples.time==BlockTimes(Block_to_analize,2)))=1;
        for i=1:length(eyedata.Events.Efix.start)
            if eyedata.Events.Efix.start(i)>BlockTimes(Block_to_analize,1) && eyedata.Events.Efix.end(i)<BlockTimes(Block_to_analize,2)
                desired_block_fixations(find(eyedata.Samples.time==eyedata.Events.Efix.start(i)):find(eyedata.Samples.time==eyedata.Events.Efix.end(i)),1)=1;
            end
        end
        

        blink_duration=mean(eyedata.Events.Eblink.duration(eyedata.Events.Eblink.duration<1000))*2.5;
        for i=1:length(eyedata.Events.Eblink.start)
            if eyedata.Events.Eblink.start(i)>BlockTimes(Block_to_analize,1) && eyedata.Events.Eblink.end(i)<BlockTimes(Block_to_analize,2) && eyedata.Events.Esacc.end(i)<BlockTimes(Block_to_analize,2) && find(eyedata.Samples.time==eyedata.Events.Eblink.start(i))>blink_duration
                desired_block_fixations(find(eyedata.Samples.time==eyedata.Events.Eblink.start(i))-blink_duration:find(eyedata.Samples.time==eyedata.Events.Eblink.end(i))+blink_duration)=0;
            end
        end
        
        saccade_duration=mean(eyedata.Events.Esacc.duration(eyedata.Events.Esacc.duration<1000))*2.5;
        for i=1:length(eyedata.Events.Esacc.start)
            if eyedata.Events.Esacc.start(i)>BlockTimes(Block_to_analize,1) && eyedata.Events.Esacc.end(i)<BlockTimes(Block_to_analize,2) && find(eyedata.Samples.time==eyedata.Events.Esacc.start(i))>saccade_duration
                desired_block_fixations(find(eyedata.Samples.time==eyedata.Events.Esacc.start(i))-saccade_duration:find(eyedata.Samples.time==eyedata.Events.Esacc.end(i))+saccade_duration)=0;
            end
        end
        
        desired_block_fixations=logical(desired_block_fixations);
        X_data=Xtemp(desired_block_fixations);
        Y_data=Ytemp(desired_block_fixations);
                
        if  strcmp(Cued_Location,'Top')
            %             Data_collection(1,(Block_to_analize-1)*120000+1:(Block_to_analize-1)*120000+length(X_data),1)=X_data;
            %             Data_collection(1,(Block_to_analize-1)*120000+1:(Block_to_analize-1)*120000+length(Y_data),2)=Y_data;
            %             Pupil_sizes(1,Block_to_analize)=mean(eyedata.Samples.pupilSize(desired_block_fixations));
            LabelsTB(1,(Block_to_analize-1)*120000+1:(Block_to_analize-1)*120000+length(X_data))=1;
            
        else
            %             Data_collection(2,(Block_to_analize-1)*120000+1:(Block_to_analize-1)*120000+length(X_data),1)=X_data;
            %             Data_collection(2,(Block_to_analize-1)*120000+1:(Block_to_analize-1)*120000+length(Y_data),2)=Y_data;
            %             Pupil_sizes(2,Block_to_analize)=mean(eyedata.Samples.pupilSize(desired_block_fixations));
            LabelsTB(1,(Block_to_analize-1)*120000+1:(Block_to_analize-1)*120000+length(X_data))=2;
            
        end
        if  strcmp(Active_Monitoring,'Active')
            Data_collection(1,(Block_to_analize-1)*120000+1:(Block_to_analize-1)*120000+length(X_data),1)=X_data;
            Data_collection(1,(Block_to_analize-1)*120000+1:(Block_to_analize-1)*120000+length(Y_data),2)=Y_data;
            Pupil_sizes(1,Block_to_analize)=mean(eyedata.Samples.pupilSize(desired_block_fixations));
            LabelsAM(1,(Block_to_analize-1)*120000+1:(Block_to_analize-1)*120000+length(X_data))=1;
        else
            Data_collection(2,(Block_to_analize-1)*120000+1:(Block_to_analize-1)*120000+length(X_data),1)=X_data;
            Data_collection(2,(Block_to_analize-1)*120000+1:(Block_to_analize-1)*120000+length(Y_data),2)=Y_data;
            Pupil_sizes(2,Block_to_analize)=mean(eyedata.Samples.pupilSize(desired_block_fixations));
            LabelsAM(1,(Block_to_analize-1)*120000+1:(Block_to_analize-1)*120000+length(X_data))=2;
        end
        [Subject Block_to_analize]       
    end 
    save(['EyedataSubj',num2str(Subject),'.mat'],'Data_collection','Pupil_sizes','LabelsAM','LabelsTB')
end
ccc
%% Plotting
clc;
clear all;
load('EyedataSubj19.mat');
figure;
plot(squeeze(Data_collection(1,:,1))+(1360./2-nanmean(nanmean(Data_collection(:,:,1)))),-(squeeze(Data_collection(1,:,2))+(768./2-nanmean(nanmean(Data_collection(:,:,2))))),'b')
hold on;
plot(squeeze(Data_collection(2,:,1))+(1360./2-nanmean(nanmean(Data_collection(:,:,1)))),-(squeeze(Data_collection(2,:,2))+(768./2-nanmean(nanmean(Data_collection(:,:,2))))),'r')

line([1360./2-15 1360./2+15],[-768./2 -768./2],'Linewidth',2,'color',[0 0 0]);
line([1360./2 1360./2],[-768./2+15 -768./2-15],'Linewidth',2,'color',[0 0 0]);

line([786-20 786+20],[-196 -196],'Linewidth',2,'color',[0 0 0]);
line([786 786],[-196-20 -196+20],'Linewidth',2,'color',[0 0 0]);
line(1360-[786-20 786+20],-(768-[196 196]),'Linewidth',2,'color',[0 0 0]);
line(1360-[786 786],-(768-[196-20 196+20]),'Linewidth',2,'color',[0 0 0]);
grid on
xlim([0 1360]);
ylim([-768 0]);
legend Active Monit
xlabel('X-position')
ylabel('Y-position')

% Plotting Histogramns
figure;
for j=1:2
    sample_data(1,:)=Data_collection(j,~isnan(Data_collection(j,:,1)),1)-nanmean(nanmean(Data_collection(:,:,1)));
    sample_data(2,:)=-(Data_collection(j,~isnan(Data_collection(j,:,2)),2)-nanmean(nanmean(Data_collection(:,:,2))));
    
    mult=50; % some constant which is used
    X_mean2=[786 -196];
    X_mean1=[1360-786 -768+196];
    fixation_cross=[1360/2 -768/2];
    total_distance=sqrt((X_mean2(1)-X_mean1(1)).^2+(X_mean2(2)-X_mean1(2)).^2);
    lent=(X_mean2-X_mean1)';
    t=[-mult mult];
    dim1=X_mean1'+t(1)*lent;
    dim2=X_mean1'+t(2)*lent;
    X_projected=zeros(size(sample_data,2),1);
    coordinates=zeros(size(sample_data,2),2);
    for i=1:size(sample_data,2)
        coordinates(i,:)=project_point_to_line_segment(dim1',dim2',sample_data(:,i)');
        if coordinates(i,1)<fixation_cross(1)
            X_projected(i,1)=-sqrt((fixation_cross(1)-coordinates(i,1)).^2+(fixation_cross(2)-coordinates(i,2)).^2);
        else
            X_projected(i,1)=sqrt((fixation_cross(1)-coordinates(i,1)).^2+(fixation_cross(2)-coordinates(i,2)).^2);
        end
    end
    if j==1
        h1=histogram((X_projected./total_distance)*14);
        data1=(X_projected./total_distance)*14;
        h1.Normalization = 'probability';
        h1.BinWidth = 0.25;
        hold on;
    else
        h2=histogram((X_projected./total_distance)*14);
        data2=(X_projected./total_distance)*14;
        h2.Normalization = 'probability';
        h2.BinWidth = 0.25;
    end
    clearvars sample_data
end
line([7 7],[0 1],'Color','blue')
line([-7 -7],[0 1],'Color','red')

if length(data1)>length(data2)
    data1=data1(1:length(data2));
else
    data2=data2(1:length(data1));
end
[h p]=ttest(data1,data2);

legend({['Active M=',num2str(mean(data1)),' S=',num2str(std(data1))],['Monit M=',num2str(mean(data2)),' S=',num2str(std(data2))],'Top dot center', 'Bot dot center'});
grid on
xlim([-10 10])
ylim([0 0.3])
txt = ['Test result p=',num2str(p)];
text(-2.5,0.2,txt)

xlabel('Deviation from fixation [dva]')
ylabel('Proportion of fixation time in block')
%% Pooled analysis

for subj=[1:14 16:19]
    load(['EyedataSubj',num2str(subj),'.mat']);

    for j=1:2
        tmp=Data_collection(j,~isnan(Data_collection(j,:,1)),1)-nanmean(nanmean(Data_collection(:,:,1)));
        sample_data(1,1:length(tmp))=tmp;
        tmp=-(Data_collection(j,~isnan(Data_collection(j,:,2)),2)-nanmean(nanmean(Data_collection(:,:,2))));
        sample_data(2,1:length(tmp))=tmp;
        
        mult=50; % some constant which is used
        X_mean2=[786 -196];
        X_mean1=[1360-786 -768+196];
        fixation_cross=[1360/2 -768/2];
        total_distance=sqrt((X_mean2(1)-X_mean1(1)).^2+(X_mean2(2)-X_mean1(2)).^2);
        lent=(X_mean2-X_mean1)';
        t=[-mult mult];
        dim1=X_mean1'+t(1)*lent;
        dim2=X_mean1'+t(2)*lent;
        X_projected=zeros(size(sample_data,2),1);
        coordinates=zeros(size(sample_data,2),2);
        for i=1:size(sample_data,2)
            coordinates(i,:)=project_point_to_line_segment(dim1',dim2',sample_data(:,i)');
            if coordinates(i,1)<fixation_cross(1)
                X_projected(i,1)=-sqrt((fixation_cross(1)-coordinates(i,1)).^2+(fixation_cross(2)-coordinates(i,2)).^2);
            else
                X_projected(i,1)=sqrt((fixation_cross(1)-coordinates(i,1)).^2+(fixation_cross(2)-coordinates(i,2)).^2);
            end
        end
        means(j,subj)=nanmean((X_projected./total_distance)*14);
        stds(j,subj)=nanstd((X_projected./total_distance)*14);       
        Pupils(j,subj)=mean(Pupil_sizes(j,Pupil_sizes(j,:)>0));
        clearvars -except subj j means stds Pupils Data_projected Data_collection Pupil_sizes
    end
    subj
end

save('ActMon1.mat','means','stds','Pupils')
h1=histogram(means(1,:));
h1.Normalization = 'probability';
h1.BinWidth = 0.25;
hold on;
h2=histogram(means(2,:));
h2.Normalization = 'probability';
h2.BinWidth = 0.25;
legend({'Active','Monit'});
grid on
xlabel('Values of means [dva]')
ylabel('Proportion of subjects')
[h1,p1]=signrank(means(1,:),means(2,:))
txt = ['Test result p=',num2str(h1)];
text(0,0.2,txt)

figure;
h1=histogram(stds(1,:));
h1.Normalization = 'probability';
h1.BinWidth = 0.25;
hold on;
h2=histogram(stds(2,:));
h2.Normalization = 'probability';
h2.BinWidth = 0.25;
legend({'Active','Monit'});
grid on
xlabel('Values of stds [dva]')
ylabel('Proportion of subjects')
[h2,p2]=signrank(stds(1,:),stds(2,:))
txt = ['Test result p=',num2str(h2)];
text(0,0.2,txt)

figure;
h1=histogram(Pupils(1,[1:14 16:end]));
h1.Normalization = 'probability';
h1.BinWidth = 300;
hold on;
h2=histogram(Pupils(2,[1:14 16:end]));
h2.Normalization = 'probability';
h2.BinWidth = 300;
legend({'Active','Monit'});
grid on
xlabel('Pupil area [px]')
ylabel('Proportion of subjects')
[h3,p3]=signrank(Pupils(1,[1:14 16:end]),Pupils(2,[1:14 16:end]))
txt = ['Test result p=',num2str(h3)];
text(2300,0.2,txt)

%% Deviation towards target compared across active vs monit
% doing anova not very successful
clc;
clear all;
close all;
    load('Order_of_blocks_for_Psych.mat');

    Data_projected{1,19}=[];
    Data_projected{2,19}=[];
    Data_projected{3,19}=[];
    
for subj=[1:10 11:14 16:19]
    load(['EyedataSubj',num2str(subj),'.mat']);
    for j=1:2
        tmp=Data_collection(j,~isnan(Data_collection(j,:,1)),1)-nanmean(nanmean(Data_collection(:,:,1)));
        sample_data(1,1:length(tmp))=tmp;
        tmp=-(Data_collection(j,~isnan(Data_collection(j,:,2)),2)-nanmean(nanmean(Data_collection(:,:,2))));
        sample_data(2,1:length(tmp))=tmp;
        
        mult=50; % some constant which is used
        X_mean2=[786 -196];
        X_mean1=[1360-786 -768+196];
        fixation_cross=[1360/2 -768/2];
        total_distance=sqrt((X_mean2(1)-X_mean1(1)).^2+(X_mean2(2)-X_mean1(2)).^2);
        lent=(X_mean2-X_mean1)';
        t=[-mult mult];
        dim1=X_mean1'+t(1)*lent;
        dim2=X_mean1'+t(2)*lent;
        X_projected=zeros(size(sample_data,2),1);
        coordinates=zeros(size(sample_data,2),2);
        
        for i=1:size(sample_data,2)
            coordinates(i,:)=project_point_to_line_segment(dim1',dim2',sample_data(:,i)');
            if coordinates(i,1)<fixation_cross(1)
                X_projected(i,1)=-sqrt((fixation_cross(1)-coordinates(i,1)).^2+(fixation_cross(2)-coordinates(i,2)).^2);
            else
                X_projected(i,1)=sqrt((fixation_cross(1)-coordinates(i,1)).^2+(fixation_cross(2)-coordinates(i,2)).^2);
            end
        end
        
        means(j,subj)=nanmean((X_projected./total_distance)*14);
        stds(j,subj)=nanstd((X_projected./total_distance)*14);
        Pupils(j,subj)=mean(Pupil_sizes(j,Pupil_sizes(j,:)>0));
        
        Data_projected{1,subj}=horzcat(Data_projected{1,subj},((X_projected./total_distance)*14)');
        Data_projected{2,subj}=horzcat(Data_projected{2,subj},j*ones(1,length(X_projected)));
        clearvars -except meansAMTB subj j means stds Pupils Data_projected Data_collection Pupil_sizes LabelsAM LabelsTB meansTB
    end
    Data_projected{3,subj}=LabelsTB(1,[(~isnan(LabelsTB) & ~isnan(Data_collection(1,:,1))) | (~isnan(LabelsTB) & ~isnan(Data_collection(2,:,1)))]);

    for AM=1:2
        for TB=1:2
            meansAMTB(AM,TB,subj)=nanmean(Data_projected{1,subj}([Data_projected{2,subj}==AM & Data_projected{3,subj}==TB]));
        end
    end
    subj
end
meansAMTB(:,:,15)=nan;

h1=histogram(meansAMTB(1,1,:),'FaceColor',[0 0.1 1]);
h1.Normalization = 'probability';
h1.BinWidth = 0.1;
hold on;
h2=histogram(meansAMTB(2,1,:),'FaceColor',[0 0.1 0.6]);
h2.Normalization = 'probability';
h2.BinWidth = 0.1;

h1=histogram(meansAMTB(1,2,:),'FaceColor',[1 0.1 0]);
h1.Normalization = 'probability';
h1.BinWidth = 0.1;
hold on;
h2=histogram(meansAMTB(2,2,:),'FaceColor',[0.6 0.1 0]);
h2.Normalization = 'probability';
h2.BinWidth = 0.1;

legend({['ActiveTop mean=',num2str(nanmean(meansAMTB(1,1,:)))],['MonitTop mean=',num2str(nanmean(meansAMTB(2,1,:)))],...
    ['ActiveBot mean=',num2str(nanmean(meansAMTB(1,2,:)))],['MonitBot mean=',num2str(nanmean(meansAMTB(2,2,:)))]});
grid on
xlabel('Deviation from center [dva]')
ylabel('Proportion of subjects')

[h1,p1]=signrank(squeeze(meansAMTB(1,1,:)),squeeze(meansAMTB(2,1,:)))
txt = ['Top p=',num2str(h1)];
text(1,0.2,txt)

[h1,p1]=signrank(squeeze(meansAMTB(1,2,:)),squeeze(meansAMTB(2,2,:)))
txt = ['Bottom p=',num2str(h1)];
text(-1.5,0.2,txt)

[h1,p1]=signrank(squeeze(meansAMTB(1,1,:)),squeeze(meansAMTB(1,2,:)))
txt = ['Active p=',num2str(h1)];
text(-0.2,0.28,txt)

[h1,p1]=signrank(squeeze(meansAMTB(2,1,:)),squeeze(meansAMTB(2,2,:)))
txt = ['Monit p=',num2str(h1)];
text(-0.2,0.25,txt)

for subj=[1:10 11:14 16:19]
    y = Data_projected{1,subj}';
    g1 = Data_projected{2,subj};
    g2 = Data_projected{3,subj};
    p(:,subj) = anovan(y,{g1,g2},'model','interaction','varnames',{'g1','g2'})
end
%% Pooled of pools

for AM=1:2
    for TB=1:2
        meansAMTB(AM,TB,subj)=Data_projected{1,subj}([Data_projected{2,subj}==AM & Data_projected{3,subj}==TB]);
    end
end




%% Functions
% Block_to_analize=5;
% load('Order_of_blocks_for_Psych.mat');
%
% if Active_blocks(Subject,Block_to_analize)==1
%     Active_Monitoring='Active'
% else
%     Active_Monitoring='Monitoring'
% end
% if Top_blocks(Subject,Block_to_analize)==1
%     Cued_Location='Top'
% else
%     Cued_Location='Bottom'
% end
%
% plot(eyedata,find(eyedata.Samples.time==BlockTimes(Block_to_analize,1)), find(eyedata.Samples.time==BlockTimes(Block_to_analize,2)))
%
% % [heatMap, gaze, plotRange]=heatmap(eyedata,find(eyedata.Samples.time==BlockTimes(Block_to_analize,1)), find(eyedata.Samples.time==BlockTimes(Block_to_analize,2)));
% % image(heatMap)
