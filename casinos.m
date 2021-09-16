% Casinos task pilot
% C. Hassall
% September, 2016
% Note: This must be able to run on the Brain and Cognition ERP Lab computers
% (MATLAB 2011 or 2013)

close all; clear variables; rng('shuffle');

just_testing = 0; % Set to 1 in order to run in windowed mode (command window visible)

% IO Stuff for the ERP Lab
if ~just_testing
    dio = digitalio('parallel',1);
    addline(dio,0:7,'out');
    putvalue(dio,0);
end

% Define control keys
KbName('UnifyKeyNames');
ExitKey = KbName('ESCAPE');
left_kb_key	= KbName('s');
right_kb_key = KbName('k');
proceed_key = KbName('p');
show_key = KbName('f');

% Get participant info
while 1
    clc;
    p_number = input('Enter the participant number:\n','s');  % get the subject name/number
    rundate = datestr(now, 'yyyymmdd-HHMMSS');
    filename = strcat('casinos_', rundate, '_', p_number, '.txt');
    checker1 = ~exist(filename,'file');
    checker2 = isnumeric(str2double(p_number)) && ~isnan(str2double(p_number));
    if checker1 && checker2
        break;
    else
        disp('Invalid number, or filename already exists.');
        WaitSecs(1);
    end
end
sex = input('Sex (M/F): ','s');
age = input('Age: ');
handedness = input('Handedness (L/R): ','s');

% Variable for behavioural data
p_data = [];

% % % Physical display properties % % %
if just_testing % (iMac)
    viewingDistance = 800; % mm, approximately
    screenWidth = 95*2; % mm
    screenHeight = 71*2; % mm
    horizontalResolution = 800; % Pixels
    verticalResolution = 600; % Pixels
else % (Holroyd Lab)
    viewingDistance = 800; % mm, approximately
    screenWidth = 375; % mm
    screenHeight = 300; % mm
    horizontalResolution = 1280; % Pixels
    verticalResolution = 1040; % Pixels
end
horizontalPixelsPerMM = horizontalResolution/screenWidth;
verticalPixelsPerMM = verticalResolution/screenHeight;


% % % Task parameters (Text, Stimuli, and Feedback) % % %

% Number of trials and blocks
% Note that each trial should take around 4-5 seconds.
num_blocks = 3; % 3 blocks total (one for each condition)
block_names = {'50-50 ', 'mixed ', '80-20 '};
if just_testing
    num_trials = 6;
    rest_break_freq = 3;
else
    num_trials = 144; % 192 trials per block (must be a multiple of 6), around 16 minutes
    rest_break_freq = 48; % Every 64 trials
end
num_stimuli = 6;
cents_per_point = 3; % Note that we're aiming for around $5 - $10 per participant

% Text properties
normal_font_size = 16;
normal_font = 'Arial';
larger_text_size = 40;
text_colour = [255 255 255];

% Fixation properties (a plus sign, written as text)
fixation_colour = [0 0 0];
fixation_size = 60;

% Stim properties
stimDegrees = 2;
stimMMs = 2 * viewingDistance * tand(stimDegrees/2);
my_colours = [255 0 0; 128 0 0; 255 255 0; 128 128 0; 0 255 0; 0 128 0; 0 255 255; 0 128 128; 0 0 255; 0 0 128; 255 0 255; 128 0 128; 0 0 0; 255 255 255; 255 192 203; 255 69 0; 255 165 0; 173 255 47];
num_colours = length(my_colours);
my_shapes = [1 2 3 4 5 6];
my_shape_names = {'circle' 'triangle' 'diamond' 'pentagon' 'hexagon' 'octagon'};
num_shapes = length(my_shapes);
my_shape_sides = [-1 3 4 5 6 8]; % Number of sides
shape_order = randperm(num_shapes);
my_shapes = my_shapes(shape_order);
colour_order = randperm(18);
my_colours = my_colours(colour_order,:);
% Setup block colour order
% For the mixed block, colours 1-3 are the 50-50 stimulis, and colours 4-6
% are the 80-20 stimuli
block_colours = repmat(1:6,num_blocks, num_trials/num_stimuli);
for b = 1:num_blocks
    new_order = randperm(num_trials);
    block_colours(b,:) = block_colours(b,new_order);
end

% Feedback properties
fb_names = {'a cherry ', 'a lemon ', 'an apple ', 'an orange ', 'a plum '};
fbDegrees = 2;
fbMMs = 2 * viewingDistance *tand(fbDegrees/2);

% Set up SR mappings
correct_responses = nan(num_blocks,num_stimuli);
for b = 1:num_blocks
    correct_responses(b,:) = Shuffle(repmat([1 2],1,num_stimuli/2));
end
p_1 = 0.5; % outcome contingency for low-value blocks
p_2 = 0.8; % outcome contingency for high-value blocks

% Outcome probabilities for each stimulus within each block type
block_ps = [repmat(p_1,1,num_stimuli); [repmat(p_1,1,num_stimuli/2) repmat(p_2,1,num_stimuli/2)]; repmat(p_2,1,num_stimuli)];

% Choose a block order and feedback type automatically based on the participant number
p_check = str2double(p_number);
if ismember(p_check,[1 7 13])
    block_types = [1 2 3];
elseif ismember(p_check,[2 8 14])
    block_types = [1 3 2];
elseif ismember(p_check, [3 9 15])
    block_types = [2 1 3];
elseif ismember(p_check, [4 10 16])
    block_types = [2 3 1];
elseif ismember(p_check, [5 11 17])
    block_types = [3 1 2];
elseif ismember(p_check, [6 2 18])
    block_types = [3 2 1];
else
    block_types = [1 2 3];
end

% Start of experiment
try
    if ~just_testing
        % Hide the cursor and "disable" keyboard
        % Press ctl-c to bring it back
        ListenChar(2);
        HideCursor();
    end
    
    % Set up a PTB window
    background_colour = [0 0 0];
    if just_testing
        [win, rec] = Screen('OpenWindow', 0, background_colour, [0 0 horizontalResolution verticalResolution],32,2); % Windowed, for testing
    else
        [win, rec] = Screen('OpenWindow', 0, background_colour); % Full screen, for experiment
    end
    Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    % Now that a window has been opened, we can define stimuli properties
    % in terms of degrees of visual angle
    [xmid, ymid] = RectCenter(rec); % Midpoint of display
    
    % Pixel properties of stim and feedback
    stimHorizontalPixels = horizontalPixelsPerMM * stimMMs;
    stimVerticalPixels = verticalPixelsPerMM * stimMMs;
    fbHorizontalPixels = horizontalPixelsPerMM * fbMMs;
    fbVerticalPixels = verticalPixelsPerMM * fbMMs;
    feedback_rect = [xmid - fbHorizontalPixels/2 ymid - fbVerticalPixels/2 xmid + fbHorizontalPixels/2 ymid + fbVerticalPixels/2];
    
    % Make the background texture
    background = imread('background.jpg');
    background = Screen('MakeTexture', win, background);
    
    % Make the feedback (fruit) textures
    [fb1, map1, alpha1] = imread('cherry.png');
    [fb2, map2, alpha2] = imread('lemon.png');
    [fb3, map3, alpha3] = imread('apple.png');
    [fb4, map4, alpha4] = imread('orange.png');
    [fb5, map5, alpha5] = imread('plum.png');
    
    fb1(:,:,4) = alpha1;
    fb2(:,:,4) = alpha2;
    fb3(:,:,4) = alpha3;
    fb4(:,:,4) = alpha4;
    fb5(:,:,4) = alpha5;
    
    fb1t = Screen('MakeTexture', win, fb1);
    fb2t = Screen('MakeTexture', win, fb2);
    fb3t = Screen('MakeTexture', win, fb3);
    fb4t = Screen('MakeTexture', win, fb4);
    fb5t = Screen('MakeTexture', win, fb5);
    
    all_textures = [fb1t fb2t fb3t fb4t fb5t];
    
    % Randomly pick a fruit for win/loss feedback
    fb_choice = randperm(5);
    win_texture = all_textures(fb_choice(1));
    loss_texture = all_textures(fb_choice(2));
    
    % Store important info for this participant
    run_line = [num2str(p_number) ', ' datestr(now) ', ' sex, ', ' handedness ', ' num2str(age) ', ' block_names(block_types) ',' fb_names{fb_choice(1:2)}];
    dlmwrite('participant_info.txt',run_line,'delimiter','', '-append');
    
    % Task Variables
    total_reward = 0;
    best_score = 0;
    num_wins = 0;
    num_losses = 0;
    
    % Instructions (what the participant sees first)
    % instructions = 'On each trial, you''re going to see two coloured squares\nTry to keep your eyes on the cross in the middle of the display at all times\nWhen the cross changes colour, choose a square by pressing either the left or right mouse button.\nEach choice results in either a win or a loss.\nOne of the coloured squares is the target - choosing it is more likely to result in a win.\nOccasionally the colours of the squares will change and there will be a new target.\nYour goal is to win as often as possible.\n(press either button to proceed)';
    instructions1 = ['In this experiment you will play slot machines in three different casinos.\nEach slot machine is represented by a coloured shape (e.g. a red square).\nEach slot machine has two arms: a left arm and a right arm.\nWhen you see the coloured shape representing a particular slot machine, you will be asked to pull either the left arm or the right arm.\nLEFT = ''s'' key, RIGHT = ''k'' key.\nAfter your choice you will be shown the outcome, in the form of fruit:\n' fb_names{fb_choice(1)} '(WIN), indicating a successful bet or ' fb_names{fb_choice(2)} '(LOSE), indicating an unsuccessful bet.\nEvery time you see ' fb_names{fb_choice(1)} '(WIN), $0.03 will be added to your total.\nEvery time you see ' fb_names{fb_choice(2)} '(LOSE), $0.00 will be added to your total.\nFor each stimulus, pulling one of the arms may be more likely to result in a WIN compared to pulling the other arm.\nYour goal is to accumulate as much money as possible.\nEach casino has its own set of coloured shapes.\nThe words NEW CASINO will be displayed when you enter a new casino and begin a new game.\n(press any key to see more instructions)'];
    instructions2 = 'Please try to minimize all head and eye movements.\nKeep your gaze on the center of the screen at all times.\nDo not respond the moment you see the coloured shape - wait until you hear a "beep".\nYou will be given regular rest breaks - please take this opportunity to close your eyes or blink as needed.\n(press any key to do some practice trials)';
    
    Screen(win,'TextFont',normal_font);
    Screen(win,'TextSize',normal_font_size);
    DrawFormattedText(win, instructions1,'center', 'center', text_colour,[],[],[],2);
    Screen('Flip',win);
    KbPressWait();
    
    DrawFormattedText(win, instructions2,'center', 'center', text_colour,[],[],[],2);
    Screen('Flip',win);
    KbPressWait();
    
    % Do some practice trials
    num_practice_trials = 20;
    for p = 1:num_practice_trials
        if rand < 0.5
            practice_stim = 1;
            this_correct_response = 1;
        else
            practice_stim = 2;
            this_correct_response = 2;
        end
        
        Screen('DrawTexture', win, background);
        Screen('TextSize',win,fixation_size);
        DrawFormattedText(win, '+','center', 'center', fixation_colour,[],[],[],2);
        Screen('Flip',win);
        WaitSecs(1);
        
        Screen('DrawTexture', win, background)
        if practice_stim == 1
            base_rect = [0 0 stimHorizontalPixels stimVerticalPixels];
            centeredRect = CenterRectOnPointd(base_rect, xmid, ymid);
            Screen('FillRect',win,[255 0 0],centeredRect);
        else
            base_rect = [0 0 stimHorizontalPixels stimVerticalPixels];
            centeredRect = CenterRectOnPointd(base_rect, xmid, ymid);
            maxDiameter = max(base_rect) * 1.01;
            Screen('FillOval',win,[0 0 255],centeredRect,maxDiameter);
        end
        Screen('Flip',win);
        
        % Wait 500 ms, and check for early response
        responded_early = 0;
        tic;
        while toc < 1
            if KbCheck
                responded_early = 1;
            end
        end
        
        % Change the colour of the fixation cross ("go cue")
        Beeper(400,0.4,0.05); % 400 Hz sine tone for 150 ms
        
        % Wait up to 2000 ms for a response
        participant_response = 0;
        start_time = GetSecs();
        ellapsed_time = 0;
        response_time = 0;
        invalid_response = 1;
        chosen_side = 0;
        if ~responded_early
            while ellapsed_time < 2
                
                if ~participant_response
                    [keyIsDown, ~, keyCode] = KbCheck();
                    response = 0;
                    if keyIsDown
                        participant_response = 1;
                        response_time = GetSecs() - start_time;
                        
                        if keyCode(left_kb_key)
                            chosen_side = 1;
                        elseif keyCode(right_kb_key)
                            chosen_side = 2;
                        end
                        
                        % Check the response
                        if chosen_side == 1 % Left
                            invalid_response = 0;
                        elseif chosen_side == 2 % Right
                            invalid_response = 0;
                        else
                            invalid_response = 1;
                        end
                    end
                end
                ellapsed_time = GetSecs() - start_time;
            end
        end
        
        % Determine outcome
        this_trial_a_winner = -1;
        this_trial_optimal = -1;
        if ~invalid_response
            if chosen_side == this_correct_response
                this_trial_optimal = 1;
                if rand < 0.8
                    this_trial_a_winner = 1;
                else
                    this_trial_a_winner = 0;
                end
            else
                this_trial_optimal = 0;
                if rand > 0.8
                    this_trial_a_winner = 1;
                else
                    this_trial_a_winner = 0;
                end
            end
        end
        
        % Back to fixation
        Screen('DrawTexture', win, background);
        Screen('TextSize',win,fixation_size);
        DrawFormattedText(win, '+','center', 'center', fixation_colour,[],[],[],2);
        Screen('Flip',win);
        WaitSecs(0.5);
        
        % Display feedback for 1000 ms
        Screen('DrawTexture', win, background);
        Screen('TextSize',win,larger_text_size);
        if responded_early
            DrawFormattedText(win, 'TOO FAST','center', 'center', text_colour,[],[],[],2);
            Screen('Flip',win);
        elseif invalid_response
            DrawFormattedText(win, 'INVALID','center', 'center', text_colour,[],[],[],2);
            Screen('Flip',win);
        elseif this_trial_a_winner
            Screen('DrawTexture', win, win_texture, [], feedback_rect);
            Screen('Flip',win);
        else
            Screen('DrawTexture', win, loss_texture, [], feedback_rect);
            Screen('Flip',win);
        end
        WaitSecs(1);
        
        [~, ~, keyCode] = KbCheck();
        if keyCode(ExitKey)
            break;
        end
    end
    
    for b = 1:num_blocks
        
        disp(['Block: ' num2str(b)]);
        
        this_block_type = block_types(b);
        this_block_colour_order = block_colours(this_block_type,:);
        these_correct_responses = correct_responses(this_block_type,:);
        this_block_ps = block_ps(this_block_type,:);
        
        % How often participant chooses the correct square in this block
        performance = 0;
        
        Screen(win,'TextFont',normal_font);
        Screen(win,'TextSize',normal_font_size);
        DrawFormattedText(win, 'NEW CASINO\nNEW COLOURED SHAPES\n(Wait for experimenter before proceeding)', 'center', 'center', text_colour,[],[],[],2);
        Screen('Flip',win);
        [~, ~, keyCode] = KbCheck();
        while ~keyCode(proceed_key)
            [~, ~, keyCode] = KbCheck();
            
            if keyCode(show_key)
                Screen(win,'TextFont',normal_font);
                Screen(win,'TextSize',normal_font_size);
                DrawFormattedText(win, 'NEW CASINO\nNEW COLOURED SHAPES\n(Wait for experimenter before proceeding)', 'center', 'center', text_colour,[],[],[],2);
                Screen('DrawTexture', win, win_texture, [], feedback_rect + [-100 100 -100 100]);
                Screen('DrawTexture', win, loss_texture, [], feedback_rect + [100 100 100 100]);
                DrawFormattedText(win, 'WIN', xmid-100-normal_font_size/2, ymid+200, text_colour,[],[],[],2);
                DrawFormattedText(win, 'LOSE', xmid+100-normal_font_size/2, ymid+200, text_colour,[],[],[],2);
                Screen('Flip',win);
            else
                Screen(win,'TextFont',normal_font);
                Screen(win,'TextSize',normal_font_size);
                DrawFormattedText(win, 'NEW CASINO\nNEW COLOURED SHAPES\n(Wait for experimenter before proceeding)', 'center', 'center', text_colour,[],[],[],2);
                Screen('Flip',win); 
            end
        end
        
        this_casino_score = 0;
        
        % Trial loop
        for t = 1:num_trials
            
            disp(['Trial: ' num2str(t)]);
            this_trial_stimulus = this_block_colour_order(t);
            this_trial_shape = my_shapes(this_trial_stimulus);
            this_trial_colour = my_colours((b - 1)*num_stimuli+this_trial_stimulus,:);
            this_correct_response = these_correct_responses(this_trial_stimulus);
            this_trial_p = this_block_ps(this_trial_stimulus);
            
            % Determine the marker code for this trial
            marker_code = 0;
            switch this_block_type
                case 1
                    marker_code = 0;
                case 2
                    if this_trial_stimulus <= 3
                        marker_code = 10;
                    else
                        marker_code = 20;
                    end
                case 3
                    marker_code = 30;
            end
            
            disp(['Stimulus: ' num2str(this_trial_stimulus)]);
            disp(['Colour: ' num2str(this_trial_colour)]);
            disp(['Response: ' num2str(this_correct_response)]);
            disp(['P(outcome): ' num2str(this_trial_p)]);
            
            % Draw crosshairs for 500 ms
            Screen('DrawTexture', win, background);
            Screen('TextSize',win,fixation_size);
            DrawFormattedText(win, '+','center', 'center', fixation_colour,[],[],[],2);
            Screen('Flip',win);
            
            % MARK ONSET OF FIXATION CROSS HERE
            if ~just_testing
                putvalue(dio,marker_code + 1);
                WaitSecs(0.004);
                putvalue(dio,0);
            end
            
            jitter_amount = (rand() - 0.5) / 5;
            WaitSecs(0.5 + jitter_amount); % Fixation cross for 400 - 600 ms
            
            % Draw the circle for 500 ms
            Screen('DrawTexture', win, background);
            % Screen('FillOval', win , this_trial_colour, oval_rect);
            if this_trial_shape == 1
                base_rect = [0 0 stimHorizontalPixels stimVerticalPixels];
                maxDiameter = max(base_rect) * 1.01;
                centeredRect = CenterRectOnPointd(base_rect, xmid, ymid);
                Screen('FillOval', win, this_trial_colour, centeredRect, maxDiameter);
            else
                num_sides = my_shape_sides(this_trial_shape);
                anglesDeg = linspace(0, 360, num_sides + 1);
                anglesRad = anglesDeg * (pi / 180);
                radius = stimHorizontalPixels/2;
                yPosVector = sin(anglesRad) .* radius + ymid;
                xPosVector = cos(anglesRad) .* radius + xmid;
                isConvex = 1;
                Screen('FillPoly', win, this_trial_colour, [xPosVector; yPosVector]', isConvex);
            end
            Screen('Flip',win);
            
            % MARK ONSET OF STIMULI HERE
            if ~just_testing
                putvalue(dio,marker_code + 2);
                WaitSecs(0.004);
                putvalue(dio,0);
            end
            
            % Wait 500 ms, and check for early response
            responded_early = 0;
            tic;
            while toc < 1
                if KbCheck
                    responded_early = 1;
                end
            end
            
            % Change the colour of the fixation cross ("go cue")
            Beeper(400,0.4,0.05); % 400 Hz sine tone for 150 ms
            
            % MARK ONSET OF GO CUE HERE
            if ~just_testing
                putvalue(dio,marker_code + 3);
                WaitSecs(0.004);
                putvalue(dio,0);
            end
            
            % Wait up to 2000 ms for a response
            participant_response = 0;
            start_time = GetSecs();
            ellapsed_time = 0;
            response_time = 0;
            invalid_response = 1;
            chosen_side = 0;
            if ~responded_early
                while ellapsed_time < 2
                    
                    if ~participant_response
                        [keyIsDown, ~, keyCode] = KbCheck();
                        response = 0;
                        if keyIsDown
                            participant_response = 1;
                            response_time = GetSecs() - start_time;
                            
                            if keyCode(left_kb_key)
                                chosen_side = 1;
                            elseif keyCode(right_kb_key)
                                chosen_side = 2;
                            end
                            
                            % Check the response
                            if chosen_side == 1 % Left
                                
                                % MARK LEFT RESPONSE HERE
                                if ~just_testing
                                    putvalue(dio,marker_code + 4);
                                    WaitSecs(0.004);
                                    putvalue(dio,0);
                                end
                                
                                invalid_response = 0;
                            elseif chosen_side == 2 % Right
                                
                                % MARK RIGHT RESPONSE HERE
                                if ~just_testing
                                    putvalue(dio,marker_code + 5);
                                    WaitSecs(0.004);
                                    putvalue(dio,0);
                                end
                                
                                invalid_response = 0;
                            else
                                invalid_response = 1;
                            end
                        end
                    end
                    ellapsed_time = GetSecs() - start_time;
                end
            end
            
            % Determine outcome
            this_trial_a_winner = -1;
            this_trial_optimal = -1;
            if ~invalid_response
                
                if chosen_side == this_correct_response
                    this_trial_optimal = 1;
                    if rand < this_trial_p
                        this_trial_a_winner = 1;
                    else
                        this_trial_a_winner = 0;
                    end
                else
                    this_trial_optimal = 0;
                    if rand > this_trial_p
                        this_trial_a_winner = 1;
                    else
                        this_trial_a_winner = 0;
                    end
                end
            end
            
            % Back to fixation
            Screen('DrawTexture', win, background);
            Screen('TextSize',win,fixation_size);
            DrawFormattedText(win, '+','center', 'center', fixation_colour,[],[],[],2);
            Screen('Flip',win);
            
            % MARK ONSET OF FIXATION CROSS HERE
            
            % Jitter the feedback delay (400 - 600 ms)
            jitter_amount = (rand() - 0.5) / 5;
            WaitSecs(0.5 + jitter_amount);
            
            % Display feedback for 1000 ms
            Screen('DrawTexture', win, background);
            Screen('TextSize',win,larger_text_size);
            if responded_early
                DrawFormattedText(win, 'TOO FAST','center', 'center', text_colour,[],[],[],2);
                Screen('Flip',win);
                
                % MARK "TOO FAST" FEEDBACK HERE
                
            elseif invalid_response
                DrawFormattedText(win, 'INVALID','center', 'center', text_colour,[],[],[],2);
                Screen('Flip',win);
                
                % MARK "INVALID" FEEDBACK
                
            elseif this_trial_a_winner
                num_wins = num_wins + 1;
                total_reward = total_reward + 1;
                this_casino_score = this_casino_score + 1;
                Screen('DrawTexture', win, win_texture, [], feedback_rect);
                Screen('Flip',win);
                
                % MARK "WIN" FEEDBACK HERE
                if ~just_testing
                    putvalue(dio,marker_code + 6);
                    WaitSecs(0.004);
                    putvalue(dio,0);
                end
                
            else
                num_losses = num_losses + 1;
                Screen('DrawTexture', win, loss_texture, [], feedback_rect);
                Screen('Flip',win);
                
                % MARK "LOSE" FEEDBACK HERE
                if ~just_testing
                    putvalue(dio,marker_code + 7);
                    WaitSecs(0.004);
                    putvalue(dio,0);
                end
                
            end
            WaitSecs(1);
            
            this_data_line = [b t this_block_type this_trial_stimulus this_trial_p this_trial_colour this_trial_shape chosen_side responded_early invalid_response response_time*1000 this_trial_a_winner this_trial_optimal];
            dlmwrite(filename,this_data_line,'delimiter', '\t', '-append');
            
            [~, ~, keyCode] = KbCheck();
            if keyCode(ExitKey)
                break;
            end
            
            if mod(t,rest_break_freq)  == 0
                break_message = ['This is a rest break.\nCurrent total for this casino: $' num2str(this_casino_score.*cents_per_point/100,'%0.2f') '\nPress any key to continue'];
                Screen(win,'TextFont',normal_font);
                Screen(win,'TextSize',normal_font_size);
                DrawFormattedText(win, break_message, 'center', 'center', text_colour,[],[],[],2);
                Screen('DrawTexture', win, win_texture, [], feedback_rect + [-100 100 -100 100]);
                Screen('DrawTexture', win, loss_texture, [], feedback_rect + [100 100 100 100]);
                DrawFormattedText(win, 'WIN', xmid-100-normal_font_size/2, ymid+200, text_colour,[],[],[],2);
                DrawFormattedText(win, 'LOSE', xmid+100-normal_font_size/2, ymid+200, text_colour,[],[],[],2);
                Screen('Flip',win);
                KbPressWait();
            end
            
        end
        
        
        block_message = ['You are now leaving this casino.\nTotal for this casino: $' num2str(this_casino_score.*cents_per_point/100,'%0.2f') '\nPress any key to continue'];
        Screen(win,'TextFont',normal_font);
        Screen(win,'TextSize',normal_font_size);
        DrawFormattedText(win, block_message, 'center', 'center', text_colour,[],[],[],2);
        Screen('Flip',win);
        KbPressWait();
        
        [~, ~, keyCode] = KbCheck();
        if keyCode(ExitKey)
            break;
        end
    end
    
catch e
    disp(['Total reward: $' num2str(total_reward.*cents_per_point/100,'%0.2f')]);
    
    % Close the Psychtoolbox window and bring back the cursor and keyboard
    Screen('CloseAll');
    ListenChar(0);
    ShowCursor();
    rethrow(e);
end

disp(['Total reward: $' num2str(total_reward.*cents_per_point/100,'%0.2f')]);

% Close the Psychtoolbox window and bring back the cursor and keyboard
Screen('CloseAll');
ListenChar(0);
ShowCursor();