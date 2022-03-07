function TrainVigParadigmBlock(Subj)
% Screen('Preference', 'SkipSyncTests',1);
% Short MATLAB example that uses the Eyelink and Psychophysics Toolboxes
% This is the example as shown in the EyelinkToolbox article in BRMIC
% Cornelissen, Peters and Palmer 2002), but updated to also work on the
% PC version of the toolbox, and uses some new routines.
%
% Adapted after "Psychtoolbox\PsychHardware\EyelinkToolbox\EyelinkDemos\
% ShortDemos\EyelinkExample.m"
%
% HISTORY
%
% mm/dd/yy
% 07/01/08 js 	redone the structure of the experiment and added
%		integration messages to the EyeLink Data Viewer software
% 07/14/08 js 	added code to set your own EDF file name before opening
%		the experiment graphics
% 07/13/10  fwc made to work with new toolbox with callback and updated to
%               enable eye image display, added "cleanup" function,
%               reenabled try-catch


commandwindow;
%% tracker settings
if ~IsOctave
    commandwindow;
else
    more off;
end

%% 0: use eye-tracker; 1: curser as a pseudo-tracker
dummymode=0;

try
    % STEP 1
    % Added a dialog box to set your own EDF file name before opening
    % experiment graphics. Make sure the entered EDF file name is 1 to 8
    % characters in length and only numbers or letters are allowed.
    prompt = {'Enter tracker EDF file name (1 to 8 letters or numbers)'};
    dlg_title = 'Create EDF file';
    num_lines= 1;
    def     = {'DEMO'};
    answer  = inputdlg(prompt,dlg_title,num_lines,def);
    %edfFile= 'DEMO.EDF'
    edfFile = answer{1};
    fprintf('EDFFile: %s\n', edfFile );
    
    % STEP 2
    % Open a graphics window on the main screen
    % using the PsychToolbox's Screen function.
    screenNumber=max(Screen('Screens'));
    [w, wRect]=Screen('OpenWindow', screenNumber, 0,[],32,2);
    Screen(w,'BlendFunction',GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    
    % STEP 3
    % Provide Eyelink with details about the graphics environment
    % and perform some initializations. The information is returned
    % in a structure that also contains useful defaults
    % and control codes (e.g. tracker state bit and Eyelink key values).
    el=EyelinkInitDefaults(w);
    
    % STEP 4
    % Initialization of the connection with the Eyelink Gazetracker.
    % exit program if this fails.
    if ~EyelinkInit(dummymode)
        fprintf('Eyelink Init aborted.\n');
        cleanup;  % cleanup function
        return;
    end
    [v,vs]=Eyelink('GetTrackerVersion');
    fprintf('Running experiment on a ''%s'' tracker.\n', vs );
    
    % open file to record data to
    i = Eyelink('Openfile', edfFile);
    if i~=0
        fprintf('Cannot create EDF file ''%s'' ', edffilename);
        cleanup;
        %         Eyelink( 'Shutdown');
        return;
    end
    
    Eyelink('command', 'add_file_preamble_text ''Recorded by EyelinkToolbox demo-experiment''');
    [width, height]=Screen('WindowSize', screenNumber);
    
    
    % STEP 5
    % SET UP TRACKER CONFIGURATION
    % Setting the proper recording resolution, proper calibration type,
    % as well as the data file content;
    Eyelink('command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, width-1, height-1);
    Eyelink('message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, width-1, height-1);
    % set calibration type.
    Eyelink('command', 'calibration_type = HV9');
    % set parser (conservative saccade thresholds)
    Eyelink('command', 'saccade_velocity_threshold = 35');
    Eyelink('command', 'saccade_acceleration_threshold = 9500');
    
    % remote mode possible add HTARGET ( head target)
    Eyelink('command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
    Eyelink('command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,AREA,GAZERES,STATUS,INPUT,HTARGET');
    % set link data (used for gaze cursor)
    Eyelink('command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT,FIXUPDATE');
    Eyelink('command', 'link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS,INPUT,HTARGET');
    
    % allow to use the big button on the eyelink gamepad to accept the
    % calibration/drift correction target
    Eyelink('command', 'button_function 5 "accept_target_fixation"');
    
    % Tell the Eyelink to send a fixation update every 50 ms
    Eyelink('command', 'fixation_update_interval = %d', 50);
    Eyelink('command', 'fixation_update_accumulate = %d', 50);
    
    % make sure we're still connected.
    if Eyelink('IsConnected')~=1
        cleanup;
        return;
    end
    
    % STEP 6
    % Calibrate the eye tracker
    % setup the proper calibration foreground and background colors
    el.backgroundcolour = 128;
    el.foregroundcolour = 0;
    % Hide the mouse cursor;
    Screen('HideCursorHelper', w);
    EyelinkDoTrackerSetup(el);
    
    
    %% STEP 7 Main expt
    % Now starts running individual trials;
    % You can keep the rest of the code except for the implementation
    % of graphics and event monitoring
    % Each trial should have a pair of "StartRecording" and "StopRecording"
    % calls as well integration messages to the data file (message to mark
    % the time of critical events and the image/interest area/condition
    % information for the trial)
    
    load('Order_of_blocks_for_Psych.mat');
    Number_of_blocks=2;
    for Block_Num=1:Number_of_blocks
        if Active_blocks(2,Block_Num)==1
            Active_Monitoring=1;
            Task='Active';
        else
            Active_Monitoring=2;
            Task='Monitoring';
        end
        
        writeMovie=0;
        pushkey=32; %spacebar
        Baseline_time=0;
        Number_of_trials_in_blocks=32;
        refresh_rate=60;
        
        if Block_Num==1
            gray=0;
            [w, wRect]=Screen('OpenWindow',screenNumber,gray);
            if writeMovie==1
                movie = Screen('CreateMovie', w, 'MyTestMovie.mov', wRect(3), wRect(4), refresh_rate, ':CodecSettings=Videoquality=0.05 Profile=2');
            end
            Screen('TextSize',w, 18);
            HideCursor;
            KbCheck;
            WaitSecs(0.1);
            GetSecs;
            
            priorityLevel=MaxPriority(w);
            Priority(priorityLevel);
            
            str=sprintf(['...  Training phase  ...']);
            message = [str '\n\n\n\n\n\n\n\n\n\n ... STOP the dots hitting the static object  ...'...
                '\n\n\n\n ... the dots usually follow their trajectories  ...'...
                '\n\n\n\n ... but sometimes fail ...'...
                '\n\n\n\n ... ONLY press space bar IF they do not follow !!!  ...'...
                '\n\n\n\n ... you have ONLY ONE button press for each dot  ...'...
                '\n\n\n\n ... dots are in two colors: you will be told which color to protect  ...'...
                '\n\n\n\n ... please try to remian still in the chin rest during the blocks  ...'...
                '\n\n\n\n ... keep your eyes on the central object ALL THE TIME  ...'];
            DrawFormattedText(w, message,'center','center',WhiteIndex(w));
            Screen('Flip', w ,0);
            GetClicks(w);
            beepStart;
        end
        
        if First_color_blocks(2,Block_Num)==1
            target_color='RED';
            block_target_color=1;
            message = ['Block #',num2str(Block_Num),'\n\n\n\n\n ---  Respond to ',target_color, ' dots  ---\n\n\n\n\n Press any button when ready !!! \n\n\n\n\n'];
            DrawFormattedText(w, message,'center','center',[255 0 0]);
        else
            target_color='GREEN';
            block_target_color=0;
            message = ['Block #',num2str(Block_Num),'\n\n\n\n\n ---  Respond to ',target_color, ' dots  ---\n\n\n\n\n Press any button when ready !!! \n\n\n\n\n'];
            DrawFormattedText(w, message,'center','center',[0 255 0 ]);
        end
        
        Screen('Flip', w ,0);
        KbWait([],2);
        WaitSecs(0.5);
        Screen('Flip', w);
        
        obstacle_radius=20;
        boundary_radius = 100;
        x_boundary=330;
        obstacle_xposition1=wRect(3)./2;
        obstacle_yposition1=wRect(4)./2;
        
        Obstacle_1 = [obstacle_xposition1-obstacle_radius, obstacle_yposition1-obstacle_radius, obstacle_xposition1+obstacle_radius, obstacle_yposition1+obstacle_radius];
        % %         Obstacle_1_ring = [obstacle_xposition1-boundary_radius, obstacle_yposition1-boundary_radius, obstacle_xposition1+boundary_radius, obstacle_yposition1+boundary_radius];
        a_Obstacle_1=obstacle_yposition1./obstacle_xposition1;
        
        Num_moving_dots=2;
        proportion_of_events=0.5;
        
        moving_dots_radius=10;
        hitting_border_distance=boundary_radius+moving_dots_radius;
        hitting_obstacle_distance=obstacle_radius+moving_dots_radius;
        
        for y=1:obstacle_yposition1
            x=round(y/a_Obstacle_1);
            if round(sqrt(x.^2+y.^2))>hitting_obstacle_distance
                x_vertical_up_1=obstacle_xposition1+y;
                y_vertical_up_1=obstacle_yposition1-x;
                x_hit_boundary_1=obstacle_xposition1-x;
                break;
            end
        end
        a_vertical_up_1=y_vertical_up_1./x_vertical_up_1;
        x_vertical_down_1=obstacle_xposition1-(x_vertical_up_1-obstacle_xposition1);
        y_vertical_down_1=obstacle_yposition1-(y_vertical_up_1-obstacle_yposition1);
        a_vertical_down_1=y_vertical_down_1./x_vertical_down_1;
        
        trial_time=200;
        refractory_period=20;
        
        decoloring=20;
        
        defl_direction=zeros(1,Num_moving_dots*Number_of_trials_in_blocks);
        defl_direction(randsample([1:length(defl_direction)],length(defl_direction)./2))=1;
        defl_direction2=zeros(1,Num_moving_dots*Number_of_trials_in_blocks);
        defl_direction2(randsample([1:length(defl_direction2)],length(defl_direction2)./2))=1;
        
        
        bias_to_ensure=(a_vertical_down_1-a_vertical_up_1).*0.1;
        
        trial_of_targets=zeros(1,Number_of_trials_in_blocks);
        if Active_Monitoring==1
            inds=randsample(Number_of_trials_in_blocks,ceil(Number_of_trials_in_blocks*0.5));
        elseif Active_Monitoring==2
            inds=randsample(Number_of_trials_in_blocks,ceil(Number_of_trials_in_blocks*0.05));
        end
        trial_of_targets(1,inds)=1;
        
        top_targets=nan*ones(round(proportion_of_events.*Num_moving_dots),Number_of_trials_in_blocks);
        for trial=1:Number_of_trials_in_blocks
            % generating indices of all event dots
            a=zeros(1,Num_moving_dots);
            index_of_event_dots=randperm(Num_moving_dots,round(proportion_of_events.*Num_moving_dots));
            
            % generating top event dots
            number_of_events=round(length(index_of_event_dots));
            if mod(trial,2)
                a(index_of_event_dots(1:number_of_events))=randsample([a_vertical_up_1+bias_to_ensure:0.006:a_vertical_up_1+2*bias_to_ensure],round(length(index_of_event_dots)));
            else
                a(index_of_event_dots(1:number_of_events))=randsample([a_vertical_down_1-2*bias_to_ensure:0.006:a_vertical_down_1-bias_to_ensure],round(length(index_of_event_dots)));
            end
            top_events(:,trial)=index_of_event_dots;
            
            if trial_of_targets(1,trial)>0
                tmp=index_of_event_dots(1:number_of_events);
                index_of_target_dots=tmp(randperm(length(tmp),round(length(tmp).*trial_of_targets(1,trial))));
                top_targets(:,trial)=index_of_target_dots;
            else
                top_targets(:,trial)=0;
            end
            
            % generating non-target dots
            tmp_a=[0.2494*10000:fix((a_vertical_up_1-bias_to_ensure)*10000) fix((a_vertical_down_1+bias_to_ensure)*10000):(1*10000)];
            ind=(find(a==0));
            a(ind(1:round(length(ind))))=randsample(tmp_a,round((1-proportion_of_events).*Num_moving_dots))./10000;               % y=ax+b
            a_original(:,trial)=a;
        end
        a=a_original;
        
        
        trial_of_targets2=zeros(1,Number_of_trials_in_blocks);
        if Active_Monitoring==1
            inds=randsample(Number_of_trials_in_blocks,ceil(Number_of_trials_in_blocks*0.5));
        elseif Active_Monitoring==2
            inds=randsample(Number_of_trials_in_blocks,ceil(Number_of_trials_in_blocks*0.05));
        end
        trial_of_targets2(1,inds)=1;
        
        top_targets2=nan*ones(round(proportion_of_events.*Num_moving_dots),Number_of_trials_in_blocks);
        for trial=1:Number_of_trials_in_blocks
            % generating indices of all event dots
            a2=zeros(1,Num_moving_dots);
            index_of_event_dots=randperm(Num_moving_dots,round(proportion_of_events.*Num_moving_dots));
            
            % generating top event dots
            number_of_events=round(length(index_of_event_dots));
            if mod(trial,2)
                a2(index_of_event_dots(1:number_of_events))=randsample([a_vertical_up_1+bias_to_ensure:0.006:a_vertical_up_1+2*bias_to_ensure],round(length(index_of_event_dots)));
            else
                a2(index_of_event_dots(1:number_of_events))=randsample([a_vertical_down_1-2*bias_to_ensure:0.006:a_vertical_down_1-bias_to_ensure],round(length(index_of_event_dots)));
            end
            top_events2(:,trial)=index_of_event_dots;
            
            if trial_of_targets2(1,trial)>0
                tmp=index_of_event_dots(1:number_of_events);
                index_of_target_dots=tmp(randperm(length(tmp),round(length(tmp).*trial_of_targets2(1,trial))));
                top_targets2(:,trial)=index_of_target_dots(1:length(index_of_target_dots));
            else
                top_targets2(:,trial)=0;
            end
            
            % generating non-target dots
            tmp_a=[0.2494*10000:fix((a_vertical_up_1-bias_to_ensure)*10000) fix((a_vertical_down_1+bias_to_ensure)*10000):(1*10000)];
            ind=(find(a2==0));
            a2(ind(1:round(length(ind))))=randsample(tmp_a,round((1-proportion_of_events).*Num_moving_dots))./10000;               % y=ax+b
            a_original2(:,trial)=a2;
        end
        a2=a_original2;
        
        
        inds_targ=zeros(Num_moving_dots,Number_of_trials_in_blocks);
        dot_color=nan*ones(Num_moving_dots,Number_of_trials_in_blocks);
        for i=1:Number_of_trials_in_blocks
            if top_targets(i)>0
                inds_targ(top_targets(i),i)=1;
            end
        end
        inds_targB=inds_targ;
        first_color_inds=randsample(find(inds_targ),sum(sum(inds_targ)./2));
        dot_color(first_color_inds)=1;
        inds_targ(first_color_inds)=nan;
        second_color_inds=find(inds_targ==1);
        dot_color(second_color_inds)=0;
        dot_color(randsample(find(isnan(dot_color)),nansum(nansum(isnan(dot_color)))./2))=1;
        dot_color(isnan(dot_color))=0;
        
        inds_targ2=zeros(Num_moving_dots,Number_of_trials_in_blocks);
        dot_color2=nan*ones(Num_moving_dots,Number_of_trials_in_blocks);
        for i=1:Number_of_trials_in_blocks
            if top_targets2(i)>0
                inds_targ2(top_targets2(i),i)=1;
            end
        end
        inds_targB2=inds_targ2;
        first_color_inds2=randsample(find(inds_targ2==1),sum(sum(inds_targ2)./2));
        dot_color2(first_color_inds2)=1;
        inds_targ2(first_color_inds2)=nan;
        second_color_inds2=find(inds_targ2==1);
        dot_color2(second_color_inds2)=0;
        dot_color2(randsample(find(isnan(dot_color2)),nansum(nansum(isnan(dot_color2)))./2))=1;
        dot_color2(isnan(dot_color2))=0;
        
        
        % Shuffling of trajectories to get the highest distance between
        % consecutive events
        modified_a_events_t=zeros(2*Number_of_trials_in_blocks,500);
        for i=1:Number_of_trials_in_blocks
            a_events(1,i)=a(top_events(i),i);
            a_events(2,i)=a2(top_events2(i),i);
        end
        a_events_t(:,1)=reshape(a_events,[1 2*Number_of_trials_in_blocks]);
        modified_a_events_t([1:2:end],1)=a_events_t([1:2:end],1);
        modified_a_events_t([2:2:end],1)=a_vertical_down_1-a_events_t([2:2:end],1)+a_vertical_up_1;
        for sh=1:500
            a_events_t(:,sh)=reshape(Shuffle(a_events),[1 2*Number_of_trials_in_blocks]);
            modified_a_events_t([1:2:end],sh)=a_events_t([1:2:end],sh);
            modified_a_events_t([2:2:end],sh)=a_vertical_down_1-a_events_t([2:2:end],sh)+a_vertical_up_1;
            for i=2:2*Number_of_trials_in_blocks
                a_event_distance(i,sh)=abs(modified_a_events_t(i,sh)-modified_a_events_t(i-1,sh));
            end
        end
        [~,inx]=max(sum(a_event_distance));
        c=0;
        for i=1:2:Number_of_trials_in_blocks*2
            c=c+1;
            a(top_events(c),c)=a_events_t(i,inx);
        end
        c=0;
        for i=2:2:Number_of_trials_in_blocks*2
            c=c+1;
            a2(top_events2(c),c)=a_events_t(i,inx);
        end
        a_original=a;
        a_original2=a2;
        
        
        
        for dot_num=1:Num_moving_dots*Number_of_trials_in_blocks
            for xt=1:obstacle_xposition1
                yt=round(a_original(dot_num).*xt);
                if sqrt((xt-obstacle_xposition1).^2+(yt-obstacle_yposition1).^2)<(boundary_radius+moving_dots_radius)
                    break;
                end
            end
            hitting_coordinates1(dot_num,:)=[xt yt];
            yhit(dot_num)=hitting_coordinates1(dot_num,2);
            for cnt=1:150
                if defl_direction(dot_num)==1
                    yhit(dot_num)=yhit(dot_num)+1;
                else
                    yhit(dot_num)=yhit(dot_num)-1;
                end
                xtt(dot_num)=round((yhit(dot_num)-hitting_coordinates1(dot_num,2)+(-1./a_original(dot_num))*hitting_coordinates1(dot_num,1))./(-1./a_original(dot_num)));
            end
            deflection_coordinates1(dot_num,:)=[xtt(dot_num) yhit(dot_num)];
        end
        y_boundary=a_original.*x_boundary;
        
        
        for dot_num=1:Num_moving_dots*Number_of_trials_in_blocks
            for xt=1:obstacle_xposition1
                xtmp=wRect(3)-xt;
                yt=round(a_original2(dot_num).*xtmp)+wRect(4)-round(wRect(3).*a_original2(dot_num));
                
                if sqrt((xtmp-obstacle_xposition1).^2+(yt-obstacle_yposition1).^2)<(boundary_radius+moving_dots_radius)
                    break;
                end
            end
            hitting_coordinates2(dot_num,:)=[xtmp yt];
            yhit(dot_num)=hitting_coordinates2(dot_num,2);
            for cnt=1:150
                if defl_direction2(dot_num)==1
                    yhit(dot_num)=yhit(dot_num)+1;
                else
                    yhit(dot_num)=yhit(dot_num)-1;
                end
                xtt(dot_num)=round((yhit(dot_num)-hitting_coordinates2(dot_num,2)+(-1./a_original2(dot_num))*hitting_coordinates2(dot_num,1))./(-1./a_original2(dot_num)));
            end
            deflection_coordinates2(dot_num,:)=[xtt(dot_num) yhit(dot_num)];
        end
        y_boundary2=a_original2.*(wRect(3)-x_boundary)+wRect(4)-round(wRect(3).*a_original2);
        
        
        
        speedx1=4.*ones(Num_moving_dots,Number_of_trials_in_blocks)-(a_original-0.2494);
        speedyd1=speedx1.*wRect(4)./wRect(3);
        speedx2=4.*ones(Num_moving_dots,Number_of_trials_in_blocks)-(a_original2-0.2494);
        speedyd2=speedx2.*wRect(4)./wRect(3);
        
        appearance_time_temp=zeros(2*Num_moving_dots,Number_of_trials_in_blocks);
        top_indices=[1:Num_moving_dots];
        top_indices2=[Num_moving_dots+1:2*Num_moving_dots];
        %         distance_between_dots_limit=(trial_time./10);
        distance_between_dots_limit=(trial_time./10);
        for tr=1:Number_of_trials_in_blocks
            %             distance_between_events_limit=randsample([fix(trial_time./2.4) fix(trial_time./2.1)],1);
            distance_between_events_limit=randsample([fix(trial_time./6) fix(trial_time./5)],1);
            appearance_time_temp([top_indices(top_events(:,tr)) top_indices2(top_events2(:,tr))],tr)=randsample([(tr-1)*trial_time+1:distance_between_events_limit:tr*trial_time-distance_between_events_limit],length([top_indices(top_events(:,tr)) top_indices2(top_events2(:,tr))]));
            for iter=1:500
                %                 tp(iter,:)=randsample([(tr-1)*trial_time+1+distance_between_dots_limit:distance_between_dots_limit:tr*trial_time-distance_between_dots_limit],sum(appearance_time_temp(:,tr)==0));
                tp(iter,:)=randsample([(tr-1)*trial_time+1+distance_between_dots_limit:distance_between_dots_limit:tr*trial_time-distance_between_dots_limit],sum(appearance_time_temp(:,tr)==0));
                tg=appearance_time_temp(appearance_time_temp(:,tr)~=0,tr);
                for i=1:length(tp(iter,:))
                    for j=1:length(tg)
                        tmpp(i,j)=abs(tp(iter,i)-tg(j));
                    end
                end
                dist_alterns(iter)=min(min(tmpp));
            end
            [~,indmax]=max(dist_alterns);
            appearance_time_temp(appearance_time_temp(:,tr)==0,tr)=tp(indmax,:);
        end
        
        appearance_time=appearance_time_temp(1:Num_moving_dots,:);
        appearance_time2=appearance_time_temp(Num_moving_dots+1:2*Num_moving_dots,:);
        
        
        y=zeros(1,Num_moving_dots*Number_of_trials_in_blocks);
        x0=zeros(Num_moving_dots,Number_of_trials_in_blocks);
        y0=zeros(Num_moving_dots,Number_of_trials_in_blocks);
        
        yy=zeros(1,Num_moving_dots*Number_of_trials_in_blocks);
        xx0=zeros(Num_moving_dots,Number_of_trials_in_blocks);
        yy0=zeros(Num_moving_dots,Number_of_trials_in_blocks);
        
        
        counter1=zeros(1,Num_moving_dots*Number_of_trials_in_blocks);
        counter1_temp=zeros(1,Num_moving_dots*Number_of_trials_in_blocks);
        
        counter2=zeros(1,Num_moving_dots*Number_of_trials_in_blocks);
        counter2_temp=zeros(1,Num_moving_dots*Number_of_trials_in_blocks);
        temp=zeros(1,Num_moving_dots*Number_of_trials_in_blocks);
        
        key_pressed1=zeros(Num_moving_dots,trial_time*(Number_of_trials_in_blocks+0.5),Number_of_trials_in_blocks);
        key_pressed2=zeros(Num_moving_dots,trial_time*(Number_of_trials_in_blocks+0.5),Number_of_trials_in_blocks);
        key_pressedTotal=zeros(1,trial_time*(Number_of_trials_in_blocks+0.5),Number_of_trials_in_blocks);
        dist_top=nan*ones(Num_moving_dots,Number_of_trials_in_blocks);
        dist_top2=nan*ones(Num_moving_dots,Number_of_trials_in_blocks);
        
        automatically_deflected_top=nan*ones(size(appearance_time));
        manually_deflected_top=nan*ones(size(appearance_time));
        beeped_top=nan*ones(size(appearance_time));
        
        automatically_deflected_top2=nan*ones(size(appearance_time));
        manually_deflected_top2=nan*ones(size(appearance_time));
        beeped_top2=nan*ones(size(appearance_time));
        
        FirstAppearTop1=zeros(Num_moving_dots,Number_of_trials_in_blocks);
        FirstAppearTop2=zeros(Num_moving_dots,Number_of_trials_in_blocks);
        FirstShadeTop1=zeros(Num_moving_dots,Number_of_trials_in_blocks);
        FirstShadeTop2=zeros(Num_moving_dots,Number_of_trials_in_blocks);
        
        distance_traj1=nan*ones(Num_moving_dots*Number_of_trials_in_blocks,trial_time*(Number_of_trials_in_blocks+0.5));
        distance_traj2=nan*ones(Num_moving_dots*Number_of_trials_in_blocks,trial_time*(Number_of_trials_in_blocks+0.5));
        
        shades1=zeros(Num_moving_dots,Number_of_trials_in_blocks);
        shades2=zeros(Num_moving_dots,Number_of_trials_in_blocks);
        
        dot_num_closest_top=1;
        dot_num_closest_top2=1;
        refractory_counter=0;
        
        % STEP 7.1
        % Sending a 'TRIALID' message to mark the start of a trial in Data
        % Viewer.  This is different than the start of recording message
        % START that is logged when the trial recording begins. The viewer
        % will not parse any messages, events, or samples, that exist in
        % the data file prior to this message.
        Eyelink('Message', 'BLOCKID %d', Block_Num);
        
        % This supplies the title at the bottom of the eyetracker display
        Eyelink('command', 'record_status_message "BLOCK %d/%d, Task %s, Target %s"', Block_Num, Number_of_blocks, Task, target_color);
        % Before recording, we place reference graphics on the host display
        % Must be offline to draw to EyeLink screen
        Eyelink('Command', 'set_idle_mode');
        % clear tracker display and draw box at center
        Eyelink('Command', 'clear_screen 0')
        Eyelink('command', 'draw_box %d %d %d %d 15', width/2-50, height/2-50, width/2+50, height/2+50);
        
        
        % STEP 7.2
        % Do a drift correction at the beginning of each trial
        % Performing drift correction (checking) is optional for
        % EyeLink 1000 eye trackers.
        %       % EyelinkDoDriftCorrection(el);
        
        % STEP 7.3
        % start recording eye position (preceded by a short pause so that
        % the tracker can finish the mode transition)
        % The paramerters for the 'StartRecording' call controls the
        % file_samples, file_events, link_samples, link_events availability
        Eyelink('Command', 'set_idle_mode');
        WaitSecs(0.05);
        
        %         % tracker starts recording!
        Eyelink('StartRecording');
        % record a few samples before we actually start displaying
        % otherwise you may lose a few msec of data
        WaitSecs(0.1);
        
        % get eye that's tracked
        % returns 0 (LEFT_EYE), 1 (RIGHT_EYE) or 2 (BINOCULAR) depending on what data is
        eye_used = Eyelink('EyeAvailable');
        if eye_used == 2
            eye_used = 1; % use the right-eye data
        end
        
        %% fixation-dependent proceeding
        % fixation on cross for a certain time to proceed to search display
        %         proceedTime = GetSecs + 10;
        totalFixTime=0;
        Fixation_assurance_time = 3;
        radius=100;
        
        % Draw fixation cross
        Screen('FillOval',w,[70 70 70], Obstacle_1);
        Screen('Flip', w);
        Eyelink('Message', 'FixationAppeared');

        while totalFixTime < Fixation_assurance_time
            
            % if using eye-tracker
            error=Eyelink('CheckRecording');
            if(error~=0)
                disp('Error in Recording');
                break;
            end
            % we need to loop over this a few times (30 is
            % randomly chosen) so that we do not miss any events
            % and to prevent any buffer overflow
            for j=1:30
                evtype = Eyelink('GetNextDataType');
                if evtype == el.FIXUPDATE
                    if Eyelink('isconnected') == el.connected % check if eye-tracker is connected
                        evt = Eyelink('getfloatdata', evtype);% get data
                        % only process if its the desired eye &&
                        % within ROI
                        if evt.eye == eye_used
                            if sqrt((evt.gavx-width/2)^2+(evt.gavy-height/2)^2)<radius
                                totalFixTime = totalFixTime + 3; % make it consistent with fixupdate rate (1/50ms)
                                if totalFixTime >= Fixation_assurance_time % if fixation beyond threshold
                                    disp('fixation done correctly!');
                                    break;
                                else % reset if no fixation
                                    totalFixTime = 0;
                                end
                            end
                        end
                    else
                        disp('Eyelink disconnected!');
                    end
                end
            end
        end
        
        %% Stimuli onset
        % mark zero-plot time in data file
        Eyelink('Message', 'Block_Onset');
        
        for total_counter=1:trial_time*(Number_of_trials_in_blocks+1)
            % Indicate the time of the picture onset
            % (this message can be used to create an interest period in EyeLink
            % Data Viewer)
            % display onset
            
            [~,~,KeyCode]=KbCheck;
            if KeyCode(pushkey)==1
                key_pressedTotal(1,total_counter)=1;
            end
            % Check recording status, stop display if error
            
            error=Eyelink('CheckRecording');
            if(error~=0)
                break;
            end
            refractory_counter=refractory_counter+1;
            
            Screen('FillOval',w,[255 255 255], Obstacle_1);
            
            %         Screen('FrameOval',w,[255 255 20], Obstacle_1_ring);
            
            
            for dot_num=1:Num_moving_dots*Number_of_trials_in_blocks
                [~,dot_num_closest_top_total(1,tr)]=nanmin([dist_top(:,tr);dist_top2(:,tr)]);
                if dot_num_closest_top_total(1,tr)<=Num_moving_dots
                    closest_top_is_left(1,tr)=1;
                else
                    closest_top_is_left(1,tr)=0;
                end
                
                tr=ceil(dot_num./Num_moving_dots);
                dot_in_trial=dot_num-(tr-1).*Num_moving_dots;
                
                if total_counter>appearance_time(dot_num)
                    counter1(dot_num)=counter1(dot_num)+speedx1(dot_num);
                    y(dot_num)=round(a(dot_num).*(counter1(dot_num)-x0(dot_num)))+y0(dot_num);
                    
                    if sum(dot_in_trial==top_events(:,tr))==1 && (isnan(automatically_deflected_top(dot_num)) || isnan(manually_deflected_top(dot_num)))
                        dist_top(dot_in_trial,tr)=sqrt((counter1(dot_num)-obstacle_xposition1).^2+(y(dot_num)-obstacle_yposition1).^2);
                    elseif sum(dot_in_trial==top_events(:,tr))==1 && (~isnan(automatically_deflected_top(dot_num)) || ~isnan(manually_deflected_top(dot_num)))
                        dist_top(dot_in_trial,tr)=sqrt((x1-obstacle_xposition1).^2+(counter1_temp(dot_num)-obstacle_yposition1).^2);
                    end
                    
                    dist_top(dist_top(:,tr)>600,tr)=nan;
                    if sum(~isnan(dist_top(:,tr)))
                        [~,dot_num_closest_top]=nanmin(dist_top(:,tr));
                    end
                    
                    if sum(dot_in_trial==top_events(:,tr))==1 && isnan(manually_deflected_top(dot_num_closest_top,tr)) && isnan(beeped_top(dot_num_closest_top,tr)) && KeyCode(pushkey)==1 && refractory_counter>refractory_period && closest_top_is_left(1,tr)==1 && (counter1((tr-1).*Num_moving_dots+dot_num_closest_top)-obstacle_xposition1)<0
                        key_pressed1(dot_num_closest_top,total_counter,tr)=1;
                        refractory_counter=0;
                        Eyelink('Message', 'KeyTop1 %d %d', tr, dot_num_closest_top);
                        if dist_top(dot_num_closest_top,tr)<hitting_border_distance && block_target_color==dot_color((tr-1).*Num_moving_dots+dot_num_closest_top)
                            a(dot_num_closest_top,tr)=-1./a(dot_num_closest_top,tr);
                            x0(dot_num_closest_top,tr)=counter1((tr-1).*Num_moving_dots+dot_num_closest_top);
                            y0(dot_num_closest_top,tr)=y((tr-1).*Num_moving_dots+dot_num_closest_top);
                            counter1_temp((tr-1).*Num_moving_dots+dot_num_closest_top)=y0(dot_num_closest_top,tr);
                            manually_deflected_top(dot_num_closest_top,tr)=1;
                            %                             beep;
                            beepCorrect;
                        else
                            beepIncorrect;
                            %                             beep;
                            InitializePsychSound(1);
                            beeped_top(dot_num_closest_top,tr)=1;
                        end
                    end
                    
                    if sum(dot_in_trial==top_events(:,tr)) && isnan(automatically_deflected_top(dot_in_trial,tr)) && isnan(manually_deflected_top(dot_in_trial,tr)) && sum(dot_in_trial==top_targets(:,tr))==0 && (sqrt((counter1(dot_num)-obstacle_xposition1).^2+(y(dot_num)-obstacle_yposition1).^2)<hitting_border_distance)
                        Eyelink('Message', 'DefTop1 %d %d', tr, dot_in_trial);
                        x0(dot_in_trial,tr)=counter1(dot_num);
                        y0(dot_in_trial,tr)=y(dot_num);
                        a(dot_num)=-1./a(dot_num);
                        automatically_deflected_top(dot_in_trial,tr)=1;
                        counter1_temp(dot_num)=y(dot_num);
                    end
                    
                    Moving_dot = [counter1(dot_num)-moving_dots_radius, y(dot_num)-moving_dots_radius, counter1(dot_num)+moving_dots_radius, y(dot_num)+moving_dots_radius];
                    if ~isnan(automatically_deflected_top(dot_num)) || ~isnan(manually_deflected_top(dot_num))
                        if defl_direction(dot_num)
                            counter1_temp(dot_num)=counter1_temp(dot_num)+speedyd1(dot_num);
                        else
                            counter1_temp(dot_num)=counter1_temp(dot_num)-speedyd1(dot_num);
                        end
                        x1=round((counter1_temp(dot_num)-y0(dot_num)+a(dot_num)*x0(dot_num))./a(dot_num));
                        Moving_dot = [x1-moving_dots_radius,counter1_temp(dot_num)-moving_dots_radius, x1+moving_dots_radius,counter1_temp(dot_num)+moving_dots_radius];
                    end
                    
                    if dot_color(dot_num)==1
                        current_dot_color=[255 0 0];
                    else
                        current_dot_color=[0 255 0];
                    end
                    
                    if sum(dot_in_trial==top_events(:,tr))==1 && Moving_dot(1)>x_boundary && Moving_dot(1)<wRect(3)-x_boundary
                        if FirstAppearTop1(dot_in_trial,tr)==0
                            Eyelink('Message', 'FirstAppearTop1 %d %d', tr, dot_in_trial);
                            FirstAppearTop1(dot_in_trial,tr)=total_counter;
                        end
                        if (isnan(automatically_deflected_top(dot_num)) && isnan(manually_deflected_top(dot_num))) && Moving_dot(1)<obstacle_xposition1+120
                            Screen('DrawLine',w,[120],x_boundary,y_boundary(dot_num),hitting_coordinates1(dot_num,1),hitting_coordinates1(dot_num,2));
                            Screen('DrawLine',w,[120],hitting_coordinates1(dot_num,1),hitting_coordinates1(dot_num,2),deflection_coordinates1(dot_num,1),deflection_coordinates1(dot_num,2));
                            Screen('FillOval', w, current_dot_color, Moving_dot);
                        elseif (~isnan(automatically_deflected_top(dot_num)) || ~isnan(manually_deflected_top(dot_num))) && ((Moving_dot(2)./Moving_dot(1))>0.2494) && ((Moving_dot(2)./Moving_dot(1))<1)
                            if Moving_dot(1)<obstacle_xposition1
                                Screen('DrawLine',w,[120],x_boundary,y_boundary(dot_num),hitting_coordinates1(dot_num,1),hitting_coordinates1(dot_num,2));
                                Screen('DrawLine',w,[120],hitting_coordinates1(dot_num,1),hitting_coordinates1(dot_num,2),deflection_coordinates1(dot_num,1),deflection_coordinates1(dot_num,2));
                            end
                            Screen('FillOval', w, current_dot_color, Moving_dot);
                        else
                            if FirstShadeTop1(dot_in_trial,tr)==0
                                Eyelink('Message', 'FirstShadeTop1 %d %d', tr, dot_in_trial);
                                FirstShadeTop1(dot_in_trial,tr)=total_counter;
                            end
                            shades1(dot_num)=shades1(dot_num)+decoloring;
                            Screen('FillOval', w, current_dot_color-shades1(dot_num), Moving_dot);
                        end
                    end
                    if (isnan(automatically_deflected_top(dot_num)) && isnan(manually_deflected_top(dot_num)))
                        distance_traj1(dot_num,total_counter)=sqrt((counter1(dot_num)-obstacle_xposition1).^2+(y(dot_num)-obstacle_yposition1).^2);
                    else
                        distance_traj1(dot_num,total_counter)=sqrt((x1-obstacle_xposition1).^2+(counter1_temp(dot_num)-obstacle_yposition1).^2);
                    end
                    %                     if sum(dot_in_trial==top_events(:,tr))==1 && distance_traj1(dot_num,total_counter)<800
                    %                         Screen('DrawLine',w,[255],1,1,hitting_coordinates1(dot_num,1),hitting_coordinates1(dot_num,2))
                    %                     end
                end
                
                if total_counter>appearance_time2(dot_num)
                    temp(dot_num)=temp(dot_num)+speedx2(dot_num);
                    counter2(dot_num)=wRect(3)-temp(dot_num);
                    yy(dot_num)=round(a2(dot_num).*(counter2(dot_num)-xx0(dot_num)))+yy0(dot_num)+wRect(4)-round(wRect(3).*a2(dot_num));
                    
                    if sum(dot_in_trial==top_events2(:,tr))==1 && (isnan(automatically_deflected_top2(dot_num)) || isnan(manually_deflected_top2(dot_num)))
                        dist_top2(dot_in_trial,tr)=sqrt((counter2(dot_num)-obstacle_xposition1).^2+(yy(dot_num)-obstacle_yposition1).^2);
                    elseif sum(dot_in_trial==top_events2(:,tr))==1 && (~isnan(automatically_deflected_top2(dot_num)) || ~isnan(manually_deflected_top2(dot_num)))
                        dist_top2(dot_in_trial,tr)=sqrt((xx1-obstacle_xposition1).^2+(counter2_temp(dot_num)-obstacle_yposition1).^2);
                    end
                    
                    if sum(~isnan(dist_top2(:,tr)))
                        [~,dot_num_closest_top2]=nanmin(dist_top2(:,tr));
                    end
                    
                    if sum(dot_in_trial==top_events2(:,tr))==1 && isnan(manually_deflected_top2(dot_num_closest_top2,tr)) && isnan(beeped_top2(dot_num_closest_top2,tr)) && refractory_counter>refractory_period && KeyCode(pushkey)==1 && closest_top_is_left(1,tr)==0 && (counter2((tr-1).*Num_moving_dots+dot_num_closest_top2)-(obstacle_xposition1))>0
                        key_pressed2(dot_num_closest_top2,total_counter,tr)=1;
                        refractory_counter=0;
                        Eyelink('Message', 'KeyTop2 %d %d', tr, dot_num_closest_top2);
                        if  dist_top2(dot_num_closest_top2,tr)<hitting_border_distance && block_target_color==dot_color2((tr-1).*Num_moving_dots+dot_num_closest_top2)
                            a2(dot_num_closest_top2,tr)=-1./a2(dot_num_closest_top2,tr);
                            xx0(dot_num_closest_top2,tr)=counter2((tr-1).*Num_moving_dots+dot_num_closest_top2);
                            yy0(dot_num_closest_top2,tr)=yy((tr-1).*Num_moving_dots+dot_num_closest_top2);
                            counter2_temp((tr-1).*Num_moving_dots+dot_num_closest_top2)=yy0(dot_num_closest_top2,tr);
                            manually_deflected_top2(dot_num_closest_top2,tr)=1;
                            %                             beep;
                            beepCorrect;
                        else
                            beepIncorrect;
                            %                             beep;
                            %                             beeped_top2(dot_num_closest_top2,tr)=1;
                        end
                    end
                    
                    if sum(dot_in_trial==top_events2(:,tr)) && isnan(automatically_deflected_top2(dot_in_trial,tr)) && isnan(manually_deflected_top2(dot_in_trial,tr)) && sum(dot_in_trial==top_targets2(:,tr))==0 && (sqrt((counter2(dot_num)-obstacle_xposition1).^2+(yy(dot_num)-obstacle_yposition1).^2)<hitting_border_distance)
                        Eyelink('Message', 'DefTop2 %d %d', tr, dot_in_trial);
                        yy0(dot_in_trial,tr)=yy(dot_num);
                        xx0(dot_in_trial,tr)=counter2(dot_num);
                        a2(dot_num)=-1./a2(dot_num);
                        automatically_deflected_top2(dot_in_trial,tr)=1;
                        counter2_temp(dot_num)=yy0(dot_in_trial,tr);
                    end
                    
                    Moving_dot = [counter2(dot_num)-moving_dots_radius,yy(dot_num)-moving_dots_radius, counter2(dot_num)+moving_dots_radius,yy(dot_num)+moving_dots_radius];
                    
                    if ~isnan(automatically_deflected_top2(dot_num)) || ~isnan(manually_deflected_top2(dot_num))
                        if defl_direction2(dot_num)
                            counter2_temp(dot_num)=counter2_temp(dot_num)+speedyd2(dot_num);
                        else
                            counter2_temp(dot_num)=counter2_temp(dot_num)-speedyd2(dot_num);
                        end
                        xx1=round((counter2_temp(dot_num)-yy0(dot_num)+a2(dot_num).*xx0(dot_num))./a2(dot_num));
                        Moving_dot = [xx1-moving_dots_radius,counter2_temp(dot_num)-moving_dots_radius,xx1+moving_dots_radius,counter2_temp(dot_num)+moving_dots_radius];
                    end
                    
                    if dot_color2(dot_num)==1
                        current_dot_color2=[255 0 0];
                    else
                        current_dot_color2=[0 255 0];
                    end
                    
                    if sum(dot_in_trial==top_events2(:,tr))==1 && Moving_dot(1)>x_boundary && Moving_dot(1)<wRect(3)-x_boundary
                        if FirstAppearTop2(dot_in_trial,tr)==0
                            Eyelink('Message', 'FirstAppearTop2 %d %d', tr, dot_in_trial);
                            FirstAppearTop2(dot_in_trial,tr)=total_counter;
                        end
                        if (isnan(automatically_deflected_top2(dot_num)) && isnan(manually_deflected_top2(dot_num))) && Moving_dot(1)>obstacle_xposition1-120
                            Screen('DrawLine',w,[120],wRect(3)-x_boundary,y_boundary2(dot_num),hitting_coordinates2(dot_num,1),hitting_coordinates2(dot_num,2));
                            Screen('DrawLine',w,[120],hitting_coordinates2(dot_num,1),hitting_coordinates2(dot_num,2),deflection_coordinates2(dot_num,1),deflection_coordinates2(dot_num,2));
                            Screen('FillOval', w,current_dot_color2, Moving_dot);
                        elseif (~isnan(automatically_deflected_top2(dot_num)) || ~isnan(manually_deflected_top2(dot_num))) && ((Moving_dot(2)./Moving_dot(1))>0.3) && (Moving_dot(2)./Moving_dot(1))<(0.9)
                            if Moving_dot(1)>obstacle_xposition1
                                Screen('DrawLine',w,[120],wRect(3)-x_boundary,y_boundary2(dot_num),hitting_coordinates2(dot_num,1),hitting_coordinates2(dot_num,2));
                                Screen('DrawLine',w,[120],hitting_coordinates2(dot_num,1),hitting_coordinates2(dot_num,2),deflection_coordinates2(dot_num,1),deflection_coordinates2(dot_num,2));
                            end
                            Screen('FillOval', w,current_dot_color2, Moving_dot);
                        else
                            if FirstShadeTop2(dot_in_trial,tr)==0
                                Eyelink('Message', 'FirstShadeTop2 %d %d', tr, dot_in_trial);
                                FirstShadeTop2(dot_in_trial,tr)=total_counter;
                            end
                            shades2(dot_num)=shades2(dot_num)+decoloring;
                            Screen('FillOval', w,current_dot_color2-shades2(dot_num), Moving_dot);
                        end
                    end
                    if (isnan(automatically_deflected_top2(dot_num)) && isnan(manually_deflected_top2(dot_num)))
                        distance_traj2(dot_num,total_counter)=sqrt((counter2(dot_num)-obstacle_xposition1).^2+(yy(dot_num)-obstacle_yposition1).^2);
                    else
                        distance_traj2(dot_num,total_counter)=sqrt((xx1-obstacle_xposition1).^2+(counter2_temp(dot_num)-obstacle_yposition1).^2);
                    end
                end
            end
            
            
            Screen('Flip',w);
            if writeMovie==1
                Screen('AddFrameToMovie', w, CenterRect([0 0 wRect(3) wRect(4)], Screen('Rect', screenNumber)), 'frontBuffer');
            end
        end
        Eyelink('Message', 'Block_Offset');
        
        Eyelink('Message', 'Start_of_Rest_Time');
        
        Screen('Flip',w);
        WaitSecs(Baseline_time);
        Eyelink('Message', 'End_of_Rest_Time');
        
        save(['Subj_',num2str(Subj),'_Blk_',num2str(Block_Num),'_ActMon_',num2str(Active_Monitoring),...
            'train.mat'],'manually_deflected_top','top_events','top_events2','top_targets','top_targets2',...
            'automatically_deflected_top','manually_deflected_top2','automatically_deflected_top2','speedx1','speedx2',...
            'moving_dots_radius','obstacle_radius','boundary_radius','distance_traj1','beeped_top','beeped_top2',...
            'distance_traj2','hitting_obstacle_distance','hitting_border_distance',...
            'key_pressed1','key_pressed2','appearance_time','appearance_time2',...
            'a_original','a_original2','a','a2','defl_direction','defl_direction2',...
            'Num_moving_dots','Number_of_trials_in_blocks','key_pressedTotal',...
            'dot_color','dot_color2','block_target_color','block_target_color','FirstAppearTop1','FirstAppearTop2',...
            'FirstShadeTop1','FirstShadeTop2');
        
        if writeMovie==1
            Screen('FinalizeMovie', movie);
        end
        
        % stop the recording of eye-movements for the current block.
        % recommended to put this right after Block_RESULT
        Eyelink('StopRecording');
        
        
        % STEP 7.7
        % Send out necessary integration messages for data analysis
        % Send out interest area information for the block
        % See "Protocol for EyeLink Data to Viewer Integration-> Interest
        % Area Commands" section of the EyeLink Data Viewer User Manual
        % IMPORTANT! Don't send too many messages in a very short period of
        % time or the EyeLink tracker may not be able to write them all
        % to the EDF file.
        % Consider adding a short delay every few messages.
        WaitSecs(0.001);
        Eyelink('Message', '!V IAREA ELLIPSE %d %d %d %d %d %s', 1, width/2-50, height/2-50, width/2+50, height/2+50,'center');
        Eyelink('Message', '!V IAREA RECTANGLE %d %d %d %d %d %s', 2, width/4-50, height/2-50, width/4+50, height/2+50,'left');
        Eyelink('Message', '!V IAREA RECTANGLE %d %d %d %d %d %s', 3, 3*width/4-50, height/2-50, 3*width/4+50, height/2+50,'right');
        Eyelink('Message', '!V IAREA RECTANGLE %d %d %d %d %d %s', 4, width/2-50, height/4-50, width/2+50, height/4+50,'up');
        Eyelink('Message', '!V IAREA RECTANGLE %d %d %d %d %d %s', 5, width/2-50, 3*height/4-50, width/2+50, 3*height/4+50,'down');
        
        % Send messages to report trial condition information
        % Each message may be a pair of trial condition variable and its
        % corresponding value follwing the '!V TRIAL_VAR' token message
        % See "Protocol for EyeLink Data to Viewer Integration-> Trial
        % Message Commands" section of the EyeLink Data Viewer User Manual
        WaitSecs(0.001);
        Eyelink('Message', '!V Block Number =  %d', Block_Num)
        
        % STEP 7.8
        % Sending a 'BLOCK_RESULT' message to mark the end of a trial in
        % Data Viewer. This is different than the end of recording message
        % END that is logged when the trial recording ends. The viewer will
        % not parse any messages, events, or samples that exist in the data
        % file after this message.
        Eyelink('Message', 'BLOCK_RESULT 0') ;
        % %         if Block_Num>1 && mod(Block_Num,4)==1
        % %             EyelinkDoTrackerSetup(el);
        % %         end
    end
    % STEP 8
    % End of Experiment; close the file first
    % close graphics window, close data file and shut down tracker
    
    Eyelink('Command', 'set_idle_mode');
    WaitSecs(0.5);
    Eyelink('CloseFile');
    
    % download data file
    try
        fprintf('Receiving data file ''%s''\n', edfFile );
        status=Eyelink('ReceiveFile');
        if status > 0
            fprintf('ReceiveFile status %d\n', status);
        end
        if 2==exist(edfFile, 'file')
            fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, pwd );
        end
    catch
        fprintf('Problem receiving data file ''%s''\n', edfFile );
    end
    
    % STEP 9
    % run cleanup function (close the eye tracker and window).
    cleanup;
    Eyelink('ShutDown');
    Screen('CloseAll');
    
    ShowCursor;
    fclose('all');
    Priority(0);
    ListenChar(0);
    %     PsychPortAudio('Close');
    return;
catch
    ShowCursor;
    fclose('all');
    Priority(0);
    cleanup;
    %     PsychPortAudio('Close');
end
end
function cleanup
% Shutdown Eyelink:
Eyelink('Shutdown');

% Close window:
sca;
commandwindow;
end

function beepStart
% MATLAB sound
fs=15000;
duration=0.1;
freq=35000;
values=0:1/fs:duration;
a=20*sin(2*pi* freq*values);
sound(a)

% InitializePsychSound(1);
% freq=44100;
% nrchannels=1;
% pahandle=PsychPortAudio('Open',[],1,0,freq,nrchannels);
% snddata= MakeBeep(250,0.1,freq);
% PsychPortAudio('FillBuffer',pahandle,snddata);
% PsychPortAudio('Start',pahandle,[],[]);
end

function beepIncorrect
% MATLAB sound
fs=20500 ;
duration=0.03;
freq=30000;
values=0:1/fs:duration;
a=20*sin(2*pi* freq*values);
sound(a)

% InitializePsychSound(1);
% freq=44100;
% nrchannels=1;
% pahandle=PsychPortAudio('Open',[],1,0,freq,nrchannels);
% snddata= MakeBeep(250,0.1,freq);
% PsychPortAudio('FillBuffer',pahandle,snddata);
% PsychPortAudio('Start',pahandle,[],[]);
end
function beepCorrect

fs=13000;
duration=0.03;
freq=40000;
values=0:1/fs:duration;
a=10*sin(2*pi* freq*values);
sound(a)

% InitializePsychSound(1);
% freq=44100;
% nrchannels=1;
% pahandle=PsychPortAudio('Open',[],1,0,freq,nrchannels);
% snddata= MakeBeep(13000,0.1,freq);
% PsychPortAudio('FillBuffer',pahandle,snddata);
% PsychPortAudio('Start',pahandle,[],[]);
end