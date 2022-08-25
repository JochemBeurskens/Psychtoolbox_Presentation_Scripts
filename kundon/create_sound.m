function [sounds, az_out, el_out] = create_sound(az, el, audDur, audrampDur,y)

% CREATE_SOUND creates a series of sound stimuli
%
% Use as:
%  [sounds, az_out, el_out] = create_sound(az, el, dur, rampdur)
%
% Where az/el are vectors that indicate the position
% dur = duration (in seconds)
% rampdur = duration of a Hanning shaped ramp (in seconds)
%
% The directional filtering is done using the interpolated HRIR-function
% sampled at 96 kHz. By construction the sampling rate of the stimuli is 
% therefore also constrained to be 96 kHz.
%
% adapted from Jan-Mathijs Schoffelen's scripts 

fsample = 96000;
nsmp    = round(audDur*fsample); % audDur in s
nramp   = round(audrampDur*fsample); % audrampDur in s

% Create a hanning window shaped ramp
ramps           = hanning(nramp*2);
ramp_on         = ramps(1:nramp)';
ramp_off        = ramps(nramp+1:end)';
taper           = [ramp_on ones(1,nsmp-nramp*2) ramp_off]';

y = y(:,1); % just take one of the channels to avoid spatial cues

% Load hrir data and check the sampling frequency
load(fullfile('normal_hrir_interpolated.mat'));
if fsample ~= fs
  error('sound stimulus has a different sampling frequency than the hrir data')
end

% Filter white noise with hrir

% get the azimuth index
[~, azid] = min(abs(azimuth_int-az));
az_out = azimuth_int(azid);

% get the elevation index
[~, elid] = min(abs(elevation_int-el));
el_out = azimuth_int(elid);
  
sounds(:,1) = taper.*filter(hrirL_int(:,azid,elid), 1, y); % left
sounds(:,2) = taper.*filter(hrirR_int(:,azid,elid), 1, y); % right


% sounds_left_far = create_sound(-10.*ones(1,80),zeros(1,80), 0.05, 0.005); % I took the values 0.05 and 0.005 from your scripts (audio duration and ramp duration)
% sounds_left = create_sound(-3.5.*ones(1,80),zeros(1,80), 0.05, 0.005);
% sounds_right_far = create_sound(10.*ones(1,80),zeros(1,80), 0.05, 0.005); % I took the values 0.05 and 0.005 from your scripts (audio duration and ramp duration)
% sounds_right = create_sound(3.5.*ones(1,80),zeros(1,80), 0.05, 0.005);
% 
% save('allsounds.mat', 'sounds_left_far', 'sounds_left','sounds_right', 'sounds_right_far');


