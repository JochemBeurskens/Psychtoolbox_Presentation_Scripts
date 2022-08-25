clear all;clc

[y,fs] = audioread('~/projects/Psychtoolbox_Presentation_Scripts/kundon/raw/ba-200ms_96k.wav');
Fs = 96000; % sampling rate @96000Hz, always check the sampling rate with Audacity before make sound!! 

audDur  = length(y)./Fs; 
audrampDur = 0.005;

az = [-10.5 -3.5 3.5 10.5];
el = 0;

sounds_left_far  = create_sound(az(1),el,audDur,audrampDur,y);
sounds_left      = create_sound(az(2),el,audDur,audrampDur,y);
sounds_right     = create_sound(az(3),el,audDur,audrampDur,y);
sounds_right_far = create_sound(az(4),el,audDur,audrampDur,y);

save('~/projects/Psychtoolbox_Presentation_Scripts/kundon/processed/allsounds_ba-200ms_order1.mat', 'sounds_left_far', 'sounds_left','sounds_right', 'sounds_right_far');

% sound(sounds_left_far,Fs); pause;
% sound(sounds_left,Fs); pause;
% sound(sounds_right,Fs); pause;
% sound(sounds_right_far,Fs); pause;

az = [-15 -5 5 15];
el = 0;

sounds_left_far  = create_sound(az(1),el,audDur,audrampDur,y);
sounds_left      = create_sound(az(2),el,audDur,audrampDur,y);
sounds_right     = create_sound(az(3),el,audDur,audrampDur,y);
sounds_right_far = create_sound(az(4),el,audDur,audrampDur,y);

save('~/projects/Psychtoolbox_Presentation_Scripts/kundon/processed/allsounds_ba-200ms_order2.mat', 'sounds_left_far', 'sounds_left','sounds_right', 'sounds_right_far');

% sound(sounds_left_far,Fs); pause;
% sound(sounds_left,Fs); pause;
% sound(sounds_right,Fs); pause;
% sound(sounds_right_far,Fs); pause;

















