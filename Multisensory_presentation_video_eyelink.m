function Multisensory_presentation_video_eyelink
%%  decide condition matrix (randomized)
clear all;
clc;
sca;
try

% Force GetSecs and WaitSecs into memory to avoid latency later on:
GetSecs;
WaitSecs(0.1);
Stopt=0;
% Force the random number generator to a unique state, necessary to
% make the trials independent
rng('shuffle');  
D = 0.001+0.0007;

% Initialize driver, request low-latency preinit:
% InitializePsychSound(1);
% 
% % Request latency mode 2, which used to be the best one in Kun's
% % measurement
% reqlatencyclass = 2;
% 
% % Requested output frequency, may need adaptation on some audio-hw:
% freq                 = 96000; % 96kHz, because of the precomputed stimuli.
% buffersize           = 0;     % Pointless to set this. Auto-selected to be optimal.
% suggestedLatencySecs = [];
% 
% % Open audio device for low-latency output:
% pahandle = PsychPortAudio('Open', [], [], reqlatencyclass, freq, 2, buffersize, suggestedLatencySecs);
% 
% % Tell driver about hardwares inherent latency, determined via calibration
% % once:
% latbias = 0;
% prelat  = PsychPortAudio('LatencyBias', pahandle, latbias);
% postlat = PsychPortAudio('LatencyBias', pahandle);

% For this code it is assumed that the input is a set of videos that are
% different in the important parameters (ie. synchrony etc.)
% PsychDebugWindowConfiguration;
% Screen('Preference', 'SkipSyncTests', 0);
% a set of videos is loaded, we have 2 different heads, each with two
% different movements and their corresponding sounds, and then we have 2
% different synchrony settings (asynchronous or synchronous). This gives a
% total of 8 different videos that are loaded.
% Next to that there are 2 locations for the sound and 2 for the videos
% Screen('Preference', 'SkipSyncTests', 1);
% Screen('Preference', 'ConserveVRAM', 64);
cd('/home_local/meduser/Desktop/Data/data/Users/jocbeu/Experiment/McGurk_effect/');



% moviename = 'V_Ge_A_Be_Final.mov';
%commented out for testing purposes
movieNames = (dir('*.mov'));
[v_ind,loc_v,loc_a]=BalanceFactors(14,1,1:4,1:4,1:4);%sync ,[0,1] ,loc_s ,1:2
% V_names = arrayfun(@(i) movieNames(v_ind(i)),1:length(v_ind),'UniformOutput',0);
trials=length(v_ind);
% load('last_trial.mat');
%% Subject information prompt
prompt    = {'Name','task: prac or formal','session','DummyMode (0 for eyelink, 1 if not)'};% TestMode = 1 with Bitsi boxes
name      = 'Subject Information';
numLines  = 1;
default   = {'test','prac','1','1'};
answer    = inputdlg(prompt, name, numLines, default);

subj.name          = answer{1}; %max length is 8 digits
subj.task          = answer{2};
subj.session       = str2double(answer{3});
dummymode          = str2double(answer{4});
name=subj.name;
task=subj.task;
session=subj.session;
savename=['last_trial_' name '_task_' task '_session_' num2str(session) '.mat'];
%% Initialize the screen etc.
AssertOpenGL; %check that the psychtoolbox version that is installed is working properly
i=1;
background=[128, 128, 128];
    
% PsychDebugWindowConfiguration;
screenid = max(Screen('Screens')); %obtain ID's for the connected screens
% Background color will be a set below:
[win, screenRect] = Screen('OpenWindow', screenid,  background);

%% initialise eyelink
% Initialize eyelink, exit program if fails 
el = EyelinkInitDefaults(win);

% Setup the proper calibration foreground and background colours
el.backgroundcolour        = BlackIndex(el.window);
el.msgfontcolour           = WhiteIndex(el.window);
el.imgtitlecolour          = WhiteIndex(el.window);
el.calibrationtargetcolour = WhiteIndex(el.window);

% for lower resolutions you might have to play around with these values
% a little. If you would like to draw larger targets on lower res
% settings please edit PsychEyelinkDispatchCallback.m and see comments
% in the EyelinkDrawCalibrationTarget function
el.calibrationtargetsize  = 1;
el.calibrationtargetwidth = 0.5;
el.targetbeep             = 0;
el.feedbackbeep           = 0;

% Update the changes to Eyelink
EyelinkUpdateDefaults(el);

% Check if eyelink has been initiated
if ~EyelinkInit(dummymode)
    fprintf('Eyelink Init aborted.\n');
    sca; 
    return;
end

% Setup eyelink file name 
edfFile = sprintf('%s%d', subj.name, subj.session);
fprintf('EDFFile: %s\n', edfFile);

% Open edf file to record data to 
status = Eyelink('Openfile', edfFile);
if status ~=0 
    fprintf('Cannot creat EDF files ''%s'' ', edfFile);
    Eyelink('shutdown');
    sca;
    return;
end

% Make sure we're still connected.
if Eyelink('IsConnected')~=1 && ~dummymode
    cleanup;
    return;
end


% Check the version
[v vs] = Eyelink('GetTrackerVersion');
fprintf('Running experiment on a ''%s'' tracker.\n', vs );

% Set the text and screen for eyelink
Eyelink('command', 'add_file_preamble_text ''Recorded by ATTENTION PROJECT''');
[width, height] = Screen('WindowSize', screenid);

% Setup tracker configuration, setting the proper recording resolution,
% proper calibration type, as well as the data file content;
Eyelink('command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, width-1, height-1);
Eyelink('message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, width-1, height-1);

% Set calibration type.
Eyelink('command', 'calibration_type = HV5');
Eyelink('command', 'generate_default_targets = YES');
% set parser (conservative saccade thresholds)
Eyelink('command', 'saccade_velocity_threshold = 35');
Eyelink('command', 'saccade_acceleration_threshold = 9500');

% Set EDF file contents using the file_sample_data and file_event_filter
% commands 
% Set link data using link_sample_data and link_event_filter
Eyelink('command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
Eyelink('command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');

% Add "HTARGET" to record possible target data for EyeLink Remote
Eyelink('command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,AREA,HTARGET,GAZERES,STATUS,INPUT');
Eyelink('command', 'link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,HTARGET,STATUS,INPUT');

% Allow to use the big button on the eyelink gamepad to accept the
% calibration/drift correction target
Eyelink('command', 'button_function 5 "accept_target_fixation"');


EyelinkDoTrackerSetup(el);
PsychPortAudio('Close'); % have to close PsychPortAudio before recording
 % Hide the mouse cursor;
Screen('HideCursorHelper', win);

%% 
% Initialize driver, request low-latency preinit:
InitializePsychSound(1);

% Force GetSecs and WaitSecs into memory to avoid latency later on:
GetSecs;
WaitSecs(0.1);

%  Screen('OpenWindow', screenid, 0); 
resp_duration=2.0;
ifi = Screen('GetFlipInterval', win); %obtain estimated flip interval of connected screens
resp_duration    = round((D+resp_duration./ifi))*ifi;
% Screen('TextFont',  win, 'Arial');
% Screen('TextStyle', win, 0);% 0=normal, 1=bold, 2=italic, 4=underline

% Set up alpha-blending for smooth (anti-aliased) lines
Screen('BlendFunction', win, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

% Realtime scheduling: Can be used if otherwise timing is not good enough.
Priority(MaxPriority(win));

% Define some variables that are related to the screen
distanceFromDisplay = 80; % distance from display, measured in the beh 1 lab: 59 cm
% vis_diam            = 0.3; % diameter of dots

ScreenResolution    = screenRect(3:4);
ScreenSize          = [53 30]; %Beh1 BenQ monitor

PixPerCm            = ScreenResolution/ScreenSize;
PixPerDeg           = distanceFromDisplay/57*PixPerCm;
[X,Y]    = RectCenter(screenRect); %obtain x,y integer point closest to centre of screen
distanceFromDisplay=distanceFromDisplay*PixPerCm;
% Define some variables that are related to the fixation cross
% Here we set the size of the arms of our fixation corss
fixCrossDimPix = 15;
% Now we set the coordinates (all relative to zero)
lineWidthPix = 2;
FixCross_coords = [X-lineWidthPix,Y-fixCrossDimPix,X+lineWidthPix,Y+fixCrossDimPix;X-fixCrossDimPix,Y-lineWidthPix,X+fixCrossDimPix,Y+lineWidthPix];
xCoords = [-fixCrossDimPix fixCrossDimPix 0 0];
yCoords = [0 0 -fixCrossDimPix fixCrossDimPix];
allCoords = [xCoords; yCoords];

%% Now comes a section which initializes the keyboard etc.
%     KbName('UnifyKeyNames');
%     Key.quit  = KbName('ESCAPE'); % to exit the function
%     Key.start = KbName('space'); % press space to start
%     if isLive == 0 % with keyboard
%         % Defining keys for input (VLAR or VRAL, balance cross subj)
%         
% %         % left hand
% %         Key.L1 = KbName('a');
% %         Key.L2 = KbName('s');
% %         Key.L3 = KbName('d');
% %         Key.L4 = KbName('f');
% %         
% %         % right hand
% %         Key.R1 = KbName('h');
% %         Key.R2 = KbName('j');
% %         Key.R3 = KbName('k');
% %         Key.R4 = KbName('l');
% 
%         % left hand
%         Key.L1 = KbName('q');
%         Key.L2 = KbName('w');
%         Key.L3 = KbName('e');
%         Key.L4 = KbName('r');
%         
%         % right hand
%         Key.R1 = KbName('u');
%         Key.R2 = KbName('i');
%         Key.R3 = KbName('o');
%         Key.R4 = KbName('p');
%         
%         % initialize the Bitsi
% %         btsi = Bitsi(''); % testing mode, with keyboard
%         
%     else % with Bitsi keypads, or MEG button boxes
%         % left hand  -- button down
%         Key.L1 = 104; % pinky
%         Key.L2 = 103; % ring
%         Key.L3 = 102; % middle
%         Key.L4 = 101; % index
%         
%         % right hand
%         Key.R1 = 97; % index finger
%         Key.R2 = 98; % middle finger
%         Key.R3 = 99; % ring finger
%         Key.R4 = 100; % pinky finger
%         
% %         btsi = Bitsi('/dev/ttyS0'); % serial port of Bitsiin beh 1 lab, MEG probably needs another address
%     end
%     RestrictKeysForKbCheck([Key.quit,Key.start,Key.L1,Key.L2,Key.L3,Key.L4,Key.R1,Key.R2,Key.R3,Key.R4]);
%     
%     % Matrix for accuracy evaluation, these are numeric matrices,
%     % the exact values of which depend - but match - the local context
%     % (i.e. isLive or ~isLive, keypad versus keyboard)
%     if subj.hands == 0 % Vis L Aud R
%         ResponseMatrix = [Key.L1 Key.L2 Key.L3 Key.L4;
%             Key.R1 Key.R2 Key.R3 Key.R4];
%     elseif subj.hands == 1
%         ResponseMatrix = [Key.R1 Key.R2 Key.R3 Key.R4;
%             Key.L1 Key.L2 Key.L3 Key.L4];
%     end
KbName('UnifyKeyNames');

esc=KbName('ESCAPE');

space=KbName('SPACE');

far_left=KbName('Z');
left=KbName('X');
right=KbName('C');
far_right=KbName('V');

RestrictKeysForKbCheck([esc,space,far_left,left,right,far_right]);
HideCursor;  
ListenChar(2);
HideCursor;
%% setting the timing parameters and giving a number of trials per block
if subj.session == 1
    startTri = 1;
%     load(sprintf('%s_%s_design.mat',subj.name,subj.task));
%     totTri   = size(design,1);
    toTri = startTri+5;
    logmatrix = zeros(8,trials);
else
    name=subj.name;
    task=subj.task;
    session=subj.session-1;
    savename_L=['last_trial_' name '_task_' task '_session_' num2str(session) '.mat'];
    load(savename_L);
%     load(sprintf('%s_%s_design.mat',subj.name,subj.task));
%     load(sprintf('%s_%03d_%s_log',subj.name,subj.session-1,subj.task), 'endTri');
    startTri = i; %no need to add 1, did this when saving

    % this prunes the pre-constructed design to contain the trials that
    % have not yet been presented, to ensure an overall balanced design
    % if the experiment is split across more than a single session.
%     totTri   = size(design,1);
    toTri = startTri+5;
end

% fixed timing parameters
precue_duration  = 0.1; % 100ms
precue_duration  = round((D+precue_duration)./ifi)*ifi; % integer number of refreshes
cue_duration     = 2.5; % 2.5s attention cue 
cue_duration     = round((D+cue_duration)./ifi)*ifi; % integer number of refreshes
fix_duration     = 1;%+0.5*rand(totTri,1); % ranging from 3~3.5s, uniform random
fix_duration     = round((D+fix_duration)./ifi)*ifi; % integer number of refreshes
postcue_duration = 0.7; % present a 700ms fixation cross aft stimulus onset
postcue_duration = round((D+postcue_duration)./ifi)*ifi; % integer number of refreshes
resp_duration    = 2.0;
resp_duration    = round((D+resp_duration)./ifi)*ifi; % integer number of refreshes

nperblock = 20; % consistent with the design matrix 

PsychImaging('PrepareConfiguration');
%% stop depends on movie
% try
    % Before recording, we place reference graphics on the host display
    % Must be offline to draw to EyeLink screen
    Eyelink('Command', 'set_idle_mode');
    % clear tracker display and draw box at center as a notification
    % for subj's performance
    Eyelink('Command', 'clear_screen 0')
    Eyelink('command', 'draw_box %d %d %d %d 15', width/2-100, height/2-100, width/2+100, height/2+100);

    % start recording eye position (preceded by a short pause so that
    % the tracker can finish the mode transition)
    % The paramerters for the 'StartRecording' call controls the
    % file_samples, file_events, link_samples, link_events availability
    Eyelink('Command', 'set_idle_mode');
    WaitSecs(0.05);
    Eyelink('StartRecording');
    % record a few samples before we actually start displaying
    % otherwise you may lose a few msec of data
    [win, screenRect] = Screen('OpenWindow', screenid,  background);
    WaitSecs(0.1);
    
    abortit = 0;
    
    % Playbackrate defaults to 1:
    rate=1;
    
    
    % Load first movie. This is a synchronous (blocking) load:
    % Return full list of movie files from directory+pattern:
    for r=1:size(movieNames,1)
        movieNames(r).name = [ pwd filesep movieNames(r).name ];
    end
     
    
    %present fixation cross
%     Screen('FillRect', win, [255,255,255], FixCross_coords', background );
%     [t_prefix, prefixonset]=Screen('Flip', win);
%     WaitSecs(5);
    % Hide the mouse cursor:
%     HideCursor;                       
    
    % Show instructions...
    tsize=30;
    Screen('TextSize', win, tsize);
    [x, y]=Screen('DrawText', win, 'Start experiment' );
%     [x, y]=Screen('DrawText', win, 'Press ESC-ape key to abort anytime.', 40, y + 10 + tsize);
%     [x, y]=Screen('DrawText', win, 'Press SPACE key when you see the discs colliding', 40, y + 10 + tsize);
%     Screen('DrawText', win, 'Press any key to start the experiment...', 40, y + 10 + tsize);
%     
%     % Flip to show the grey screen:
    Screen('Flip',win); 
%     
    % Wait for keypress + release...
    KbStrokeWait;
    
    % Show cleared screen...
    Screen('Flip',win);
    % Wait a second...
    WaitSecs(0.5);
    Moment_before_loop=GetSecs;
    % Main trial loop: Do 'trials' trials...
    for i=startTri:toTri
        Eyelink('Message', 'START_TRIALID_%d', i); % for data viewer segregation 
        m_num=v_ind(i);
        iteration = v_ind(i);
        moviename=movieNames(mod(iteration, size(movieNames,1))+1).name;   
        %now setting the location, as given by tan(theta)=opposite/adjacent
        if loc_v(i) == 1
            angle=15.0;
            v_loc='far_right';
        elseif loc_v(i) == 2
            angle=5.0;
            v_loc='right';
        elseif loc_v(i) == 3
            angle=-5.0;
            v_loc='left';
        elseif loc_v(i) == 4
            angle=-15.0;
            v_loc='far_left';
        end        
        if loc_a(i) == 1
            a_loc='far_right';
        elseif loc_a(i) == 2
            a_loc='right';
        elseif loc_a(i) == 3
            a_loc='left';
        elseif loc_a(i) == 4
            a_loc='far_left';
        end   
        sync=1;
        combi_a_v=1;%these two should be coded numebers
        angle=pi*(angle/180);
        loccntr_target=tan(angle)*distanceFromDisplay;
        SizeTarget=150;
        Coords_cntr = [X+(loccntr_target-SizeTarget) Y-SizeTarget X+(loccntr_target+SizeTarget) Y+SizeTarget];
        logmatrix(1,i)=m_num;
        logmatrix(2,i)=loc_v(i);
        logmatrix(3,i)=loc_a(i);
        logmatrix(4,i)=sync;
        logmatrix(5,i)=combi_a_v;
        % Open the moviefile and query some infos like duration, framerate,
        % width and height of video frames. We could also query the total count of frames in
        % the movie, but computing 'framecount' takes long, so avoid to query
        % this property if you don't need it!
        %present fixation cross
        Screen('FillRect', win, [255,255,255], FixCross_coords', background );
        [t_prefix prefixonset prefix_flip]=Screen('Flip', win); %this flip is used as a timing marker
        Eyelink('Message', 'pre-cue fixation onset');
%         disp(t_prefix-Moment_before_loop)
%         disp(prefixonset-Moment_before_loop)
%         disp(prefix_flip-Moment_before_loop)
        [movie movieduration fps] = Screen('OpenMovie', win, moviename);
%         movieduration     = round(movieduration./ifi)*ifi;
        timeOfEvent=prefixonset+fix_duration; %this should be the time of the start of the video
        timeOfEvent=round(D+timeOfEvent./ifi)*ifi; %make it an integer number of screen flips
        startResponse=timeOfEvent+movieduration;
        startResponse=round(D+startResponse./ifi)*ifi; %make it an integer number of screen flips
        % We estimate framecount instead of querying it - faster:
        framecount = movieduration * fps;
        
        % Start playback of the movie:
        % Play 'movie', at a playbackrate = 1 (normal speed forward),
        % play it once, aka with loopflag = 0,
        % play audio track at volume 1.0  = 100% audio volume.
        Screen('PlayMovie', movie, 1, 0, 1.0);
        
        % Video playback and key response RT collection loop:
        % This loop repeats until either the subject responded with a
        % keypress to indicate s(he) detected the event in the vido, or
        % until the end of the movie is reached.
        movietexture=0;     % Texture handle for the current movie frame.
        reactiontime=-1;    % Variable to store reaction time.
        lastpts=0;          % Presentation timestamp of last frame.
        onsettime=-1;       % Realtime at which the event was shown to the subject.
        rejecttrial=0;      % Flag which is set to 1 to reject an invalid trial.
        
        WaitSecs(fix_duration);        
        frames_passed=0;
        while(movietexture>=0 && reactiontime==-1)
            % Check if a new movie video frame is ready for visual
            % presentation: This call polls for arrival of a new frame. If
            % a new frame is ready, it converts the video frame into a
            % Psychtoolbox texture image and returns a handle in
            % 'movietexture'. 'pts' contains a so called presentation
            % timestamp. That is the time (in seconds since start of movie)
            % at which this video frame should be shown on the screen.
            % Arrival of textures is automatically synchronized to the
            % audio track and to real-world time. If the video display loop
            % can't keep up with the flow of time and the soundtrack,
            % the engine will automatically skip/drop frames to keep video
            % in sync with audio as good as possible. If the pts of a new
            % texture is greater than the 'timeOfEvent' then you'll know
            % that this texture will show the visual target event as soon
            % as you draw and 'Flip' it.
            % In case that no new video texture is available yet for
            % presentation, this function will return a zero texture handle
            % to indicate this. If no new texture will become available
            % anymore, because the end of the movie is reached, it will
            % return a handle of -1 to indicate end of playback.
            
            % The 0 - flag means: Don't wait for arrival of new frame, just
            % return a zero or -1 'movietexture' if none is ready.
            
            [movietexture, pts] = Screen('GetMovieImage', win, movie, 0);
            start=GetSecs;
            
%             disp('start')
%             disp(start)
            % Is it a valid texture?
            if (movietexture>0)
                % Yes. Draw the texture into backbuffer:
                Screen('DrawTexture', win, movietexture, [], Coords_cntr,[],[],1);
                Screen('FillRect', win, [255,255,255], FixCross_coords', background );
%                 Screen('DrawLines', win, allCoords, lineWidthPix, [128 128 128]);
               
                % Flip the display to show the image at next retrace:
                % vbl will contain the exact system time of image onset on
                % screen: This should be accurate in the sub-millisecond
                % range.
                time_present=start+((timeOfEvent-Moment_before_loop+pts));
%                 disp('a');
%                 disp(start-GetSecs);
%                 disp(timeOfEvent+pts-GetSecs);
%                 s=GetSecs;
%                 disp(start-s);
%                 Eyelink('Message', 'STIMULUS_FRAME_%d',pts);
                vbl=Screen('Flip', win);%,time_present); %cannot change this timing, as this screen function in combination with the above automatically plays the video in synchronous order
                Eyelink('Message', 'Stimulus_presentation');
                %                 disp('b');
%                 disp(timeOfEvent-vbl);
                % Is this the first frame in the video?
                if (onsettime==-1 && pts >= (timeOfEvent-Moment_before_loop))
                    % Yes: This is the first frame with a pts timestamp that is
                    % equal or greater than the timeOfEvent, so 'vbl' is
                    % the exact time when the event was presented to the
                    % subject. Define it as onsettime:
                    onsettime = vbl;
                    
                    % Compare current pts to last one to see if the movie
                    % decoder skipped a frame at this crucial point in
                    % time. That would invalidate this trial.
                    if (pts - lastpts > 1.5*(1/fps))
                        % Difference to last frame is more than 1.5 times
                        % the expected difference under assumption 'no
                        % skip'. We skipped in the wrong moment!
                        rejecttrial=1;
                    end;
                end;
                 
                % Keep track of the frames pts in order to check for skipped frames:
                lastpts=pts;
                
                % Delete the texture. We don't need it anymore:
                Screen('Close', movietexture);
                movietexture=0;
            end;
               
            % Done with drawing. Check the keyboard for subjects response:
            [keyIsDown, secs, keyCode]=KbCheck;
            if (keyIsDown==1)
                % Abort requested?
                if keyCode(esc)
                    rejecttrial=5; %last trial should be rejected due to participant aborting program
                    i=i+1; %now the trial where the next session should start should be i+1, so that not twice the same trial is performed
                    logmatrix(6,i)=rejecttrial;
                    name=subj.name;
                    save(savename,'i','name','logmatrix')
                    % This signals abortion:
                    rejecttrial=-1;
                    Screen('CloseMovie', movie);
                    WaitSecs(0.8);
                    
                    % Break out of display loop:
                    eyelink_close(edfFile);
                    ShowCursor;
                    ple
                    sca
                    Priority(0);
                    PsychPortAudio('Close');
                    RestrictKeysForKbCheck([]);
                    ListenChar(0);
                    sca;  
                    psychrethrow(psychlasterror);
                    %continue; %use this in the case where we want to have
                    %the videos stopped at the moment of keypress by the
                    %subject
                else
                    % Reject this trial:
                    rejecttrial=2;
                end;                
            end;
            frames_passed=frames_passed+1;
        end; % ...of display loop...
        sec_prvs=0;
        %show the fixation cross, color coded to the response modality, for
        %resp_duration seconds
        ini_resp = GetSecs;
        Stopt = ini_resp + resp_duration;%timeOfEvent+resp_duration%GetSecs + resp_duration;
        Screen('FillRect', win, [128,218,128 ], FixCross_coords', background );
        [t_postfix, postfixonset]=Screen('Flip', win);
        Eyelink('Message', 'RESPONSE_ONSET');
        while GetSecs<Stopt
%             [realWakeupTimeSecs] = WaitSecs(2);
            % Check the state of the keyboard.
%             [ keyIsDown, seconds, keyCode ] = KbCheck;
            
%             [keyIsDown, secs, keyCode]=KbCheck;
%             WaitSecs(2);
            [keyIsDown, seconds, keyCode, deltarateSecs] = KbCheck;
            % If the user is pressing a key, then display its code number and name.
            if keyIsDown%keyCode(space) || keyCode(far_left) || keyCode(left) || keyCode(right) || keyCode(far_right)
                if keyCode(far_left) 
                    logmatrix(7,i)=1;
                end
                if keyCode(left) 
                    logmatrix(7,i)=2;
                end
                if keyCode(right) 
                    logmatrix(7,i)=3;
                end
                if keyCode(far_right) 
                    logmatrix(7,i)=4;
                end
                timeEndVideo=vbl+movieduration;
                reactiontime=seconds-startResponse; %compare the moment of keypress with the moment the video ended/last frame was shown
                logmatrix(8,i)=reactiontime;
                % Response too early (before event happened?)
                if (reactiontime<0)
                    % Reject this trial:
                    rejecttrial=2;
%                 else
%                     Valid response: Difference between 'secs' and
%                     'onsettime' is the reaction time: 
                end;
            end;
            if keyCode(esc) %|| seconds<prefixonset+movieduration

                % Note that we use find(keyCode) because keyCode is an array.
                % See 'help KbCheck'
                fprintf('You pressed key %i which is %s\n', find(keyCode), KbName(keyCode));

                if keyCode(esc)
                    rejecttrial=5; %last trial should be rejected due to participant aborting program
                    i=i+1; %now the trial where the next session should start should be i+1, so that not twice the same trial is performed
                    save(savename,'i','name','logmatrix')

                    rejecttrial=-1;
                    Screen('CloseMovie', movie);
                    WaitSecs(0.8);
                    
                    % Break out of display loop:
                    eyelink_close(edfFile);
                    ShowCursor;
                    ple
                    sca
                    Priority(0);
                    PsychPortAudio('Close');
                    RestrictKeysForKbCheck([]);
                    ListenChar(0);
                    sca;  
                    psychrethrow(psychlasterror);
                    break;
                end

                % If the user holds down a key, KbCheck will report multiple events.
                % To condense multiple 'keyDown' events into a single event, we wait until all
                % keys have been released.
                KbReleaseWait;
            end
        end 
            
        
        % Stop movie playback, in case it isn't already stopped. We do this
        % by selection of a playback rate of zero: This will also return
        % the number of frames that had to be dropped to keep audio, video
        % and realtime in sync.
        droppedcount = Screen('PlayMovie', movie, 0, 0, 0);
        if (droppedcount > 1)%0.2*framecount)
            % any frames skipped?!? Playback problems! We
            % reject this trial...
            % Psychtoolbox works with 20% as the cap
            rejecttrial=4;
        end;
        
        % Close the moviefile.
        Screen('CloseMovie', movie);

        %now present in the case of practice, whether the correct answer
        %was given
        if subj.task=='prac' 
           tsize=30;
           Screen('TextSize', win, tsize);
           txt1=['Visual target at: '  num2str(v_loc) '.' ];
           txt2=['Auditory target at: '  num2str(a_loc) '.'];
           Screen('DrawText', win, txt1 ,X-150,Y ); %note that capital X and Y are the center of the screen where the fix. cross is presented
           Screen('DrawText', win, txt2 ,X-150,Y+30 );
           %            [x, y]=Screen('DrawText', win, 'The visual target was presented at:' + str(angle) );
           Screen('Flip', win);
           WaitSecs(2.0);
           Screen('Flip', win);
        end
        
        if (reactiontime==-1 && rejecttrial==0)
            rejecttrial=3;
        end;
        
        % Print out trials result if it was a valid trial:
        logmatrix(6,i)=rejecttrial;
        if (rejecttrial==0)
            fprintf('Trial %i valid: Reaction time was %f msecs.\n', i, 1000 * reactiontime);
        end;
        
        if (rejecttrial==1)
            fprintf('Trial %i rejected due to skip in video playback at time of event.\n', i);
        end;
        
        if (rejecttrial==2)
            fprintf('Trial %i rejected. False detection by subject.\n', i);
        end;
        
        if (rejecttrial==3)
            fprintf('Trial %i rejected. No detection by subject. Asleep?!?\n', i);
        end;
        
        if (rejecttrial==4)
            fprintf('Trial %i rejected. Way too many skips in movie playback!!!\n', i);
        end;
        
        %now check if a break is needed (I set it to once every 20 trials)
        modulus=mod(i,3);
        if modulus==0
            tsize=30;
            Screen('TextSize', win, tsize);
            Screen('DrawText', win, 'Take a short break if you need it,',X-150,Y );
            Screen('DrawText', win, 'to continue press the space bar.',X-150,Y+30 );
            Screen('DrawText', win, 'To quit press the escape key.',X-150,Y+60 );
             
            % Flip to show the grey screen with wait text on it:
            Screen('Flip',win); 
            % Wait for keypress + release...
            [secs, keyCode, deltaSecs] = KbStrokeWait;
            if keyCode(esc)
            %The two lines below are not needed here, as run is aborted in waiting period    
%                 rejecttrial=5; %last trial should be rejected due to participant aborting program
%                 i=trials+1; %now the trial where the next session should start should be i+1, so that not twice the same trial is performed

                % This signals abortion:
                rejecttrial=-1;
                WaitSecs(0.8);

                %save data
                save(savename,'i','name','logmatrix')

                % Break out of display loop:
                eyelink_close(edfFile);
                ShowCursor;
                ple
                sca
                Priority(0);
                PsychPortAudio('Close');
                RestrictKeysForKbCheck([]);
                ListenChar(0);
                sca;  
                psychrethrow(psychlasterror);
                %continue; %use this in the case where we want to have
                %the videos stopped at the moment of keypress by the
                %subject
            end;       
            % Show cleared screen...
            Screen('Flip',win);
            WaitSecs(1.0);
        end
        
        % Check if aborted.
        if (rejecttrial==-1)
            % Break out of trial loop
            break;
        end;
        
        % Wait for subject to release keys:
        KbReleaseWait;
    end; % Trial done. Next trial...
    
    % Done with the experiment. Close onscreen window and finish.
    Priority(0);
    eyelink_close(edfFile);
    PsychPortAudio('Close');
    RestrictKeysForKbCheck([]);
    ListenChar(0);
    ShowCursor;  
    sca;
    fprintf('Done. Bye!\n');
    return;
catch %#ok<CTCH>
    % Error handling: Close all windows and movies, release all ressources.
%     eyelink_close(edfFile);
    ShowCursor;
    ple
    sca;
    Priority(0);
    PsychPortAudio('Close');
    RestrictKeysForKbCheck([]);
    ListenChar(0);
%     psychrethrow(psychlasterror);
end;
end