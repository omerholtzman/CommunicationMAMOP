%% Clear command window and close any figures

clc;
close all;

%% Load configuration information

PS2000aConfig;

%% Device connection

% Check if an Instrument session using the device object |ps2000aDeviceObj|
% is still open, and if so, disconnect if the User chooses 'Yes' when prompted.
if (exist('ps2000aDeviceObj', 'var') && ps2000aDeviceObj.isvalid && strcmp(ps2000aDeviceObj.status, 'open'))
    
    openDevice = questionDialog(['Device object ps2000aDeviceObj has an open connection. ' ...
        'Do you wish to close the connection and continue?'], ...
        'Device Object Connection Open');
    
    if (openDevice == PicoConstants.TRUE)
        
        % Close connection to device.
        disconnect(ps2000aDeviceObj);
        delete(ps2000aDeviceObj);
        
    else

        % Exit script if User selects 'No'.
        return;
        
    end
    
end

% Create a device object. 
% The serial number can be specified as a second input parameter.
ps2000aDeviceObj = icdevice('picotech_ps2000a_generic.mdd');

% Connect device object to hardware.
connect(ps2000aDeviceObj);

%% Obtain Signalgenerator group object
% Signal Generator properties and functions are located in the Instrument
% Driver's Signalgenerator group.

sigGenGroupObj = get(ps2000aDeviceObj, 'Signalgenerator');
sigGenGroupObj = sigGenGroupObj(1);

%% Turn off signal generator
% Sets the output to 0 V DC.

[status.setSigGenOff] = invoke(sigGenGroupObj, 'setSigGenOff');

%% Arbitrary waveform generator - set parameters
% Set parameters (2000 mVpp, 0 mV offset, 2000 Hz frequency) and define an
% arbitrary waveform.

% Configure property value(s).
plaster_freq = 1.024;
hz = 1;
set(ps2000aDeviceObj.Signalgenerator(1), 'startFrequency', 1.0 / plaster_freq * hz);
set(ps2000aDeviceObj.Signalgenerator(1), 'stopFrequency', 1.0 / plaster_freq * hz);
set(ps2000aDeviceObj.Signalgenerator(1), 'offsetVoltage', 0.0);
set(ps2000aDeviceObj.Signalgenerator(1), 'peakToPeakVoltage', 2000.0);

%% 
% Define an Arbitrary Waveform - values must be in the range -1 to +1.
% Arbitrary waveforms can also be read in from text and csv files using
% <matlab:doc('dlmread') |dlmread|> and <matlab:doc('csvread') |csvread|>
% respectively or use the |importAWGFile| function from the <https://uk.mathworks.com/matlabcentral/fileexchange/53681-picoscope-support-toolbox PicoScope
% Support Toolbox>.
%
% Any AWG files created using the PicoScope 6 application can be read using
% the above method.
bit_amount = 256;

awgBufferSize = get(sigGenGroupObj, 'awgBufferSize');
x = linspace(0, 2 * pi, awgBufferSize / bit_amount);
zero_bit = zeros(1, awgBufferSize / bit_amount);
frequency = 341.4 * 1e3;
one_bit = sin(x * frequency) * 2;  % used to be 2

imdata = imread('wasah.bmp');
bits = double(imdata);
bits = reshape(bits, bit_amount, []);
bits = bits(2:bit_amount)';
% bits = randi([0 1],1, bit_amount - 1);
y = [];
y = cat(2, y, one_bit / 3);
for bit=bits
    if bit == 0
        y = cat(2, y, zero_bit);
    end
    if bit == 1 || bit == 255
        y = cat(2, y, one_bit);
    end
end

%% Arbitrary waveform generator - simple
% Output an arbitrary waveform with constant frequency (defined above).

% Arb. Waveform : y (defined above)

[status.setSigGenArbitrarySimple] = invoke(sigGenGroupObj, 'setSigGenArbitrarySimple', y);

%%
% 
% <<../images/ps2000a_arbitrary_waveform.PNG>>
% 
% Block data acquisition properties and functions are located in the 
% Instrument Driver's Block group.

blockGroupObj = get(ps2000aDeviceObj, 'Block');
blockGroupObj = blockGroupObj(1);

% Set pre-trigger and post-trigger samples as required - the total of this
% should not exceed the value of |maxSamples| returned from the call to
% |ps2000aGetTimebase2()|. The default of 0 pre-trigger and 8192 post-trigger
% samples is used in this example.

% set(ps2000aDeviceObj, 'numPreTriggerSamples', 0);
set(ps2000aDeviceObj, 'numPostTriggerSamples', 8192 * 125 / hz);

%%
% This example uses the |runBlock()| function in order to collect a block of
% data - if other code needs to be executed while waiting for the device to
% indicate that it is ready, use the |ps2000aRunBlock()| function and poll
% the |ps2000aIsReady()| function.

% Capture a block of data:
%
% segment index: 0 (The buffer memory is not segmented in this example)

[status.runBlock] = invoke(blockGroupObj, 'runBlock', 0);

% Retrieve data values:

startIndex              = 0;
segmentIndex            = 0;
downsamplingRatio       = 1;
downsamplingRatioMode   = ps2000aEnuminfo.enPS2000ARatioMode.PS2000A_RATIO_MODE_DECIMATE;

% Provide additional output arguments for other channels e.g. chC for
% channel C if using a 4-channel PicoScope.
[numSamples, overflow, chA, chB] = invoke(blockGroupObj, 'getBlockData', startIndex, segmentIndex, ...
                                            downsamplingRatio, downsamplingRatioMode);

%% Process data
% In this example the data values returned from the device are displayed in
% plots in a Figure.

figure1 = figure('Name','PicoScope 2000 Series (A API) Example - Block Mode Capture', ...
    'NumberTitle','off');

% Calculate sampling interval (nanoseconds) and convert to milliseconds.
% Use the |timeIntervalNanoSeconds| output from the |ps2000aGetTimebase2()|
% function or calculate it using the main Programmer's Guide.
% Take into account the downsampling ratio used.
timeIntervalNanoseconds = 1;
timeNs = double(timeIntervalNanoseconds) * downsamplingRatio * double(0:numSamples - 1);
timeMs = timeNs / 1e6;

% encode bits:
slice_per_bit = 10;
% slice the B channel:
slice_size = length(chB) / (bit_amount * slice_per_bit);
number_of_slices = length(chB) / slice_size;
bit_slices = ones(1, number_of_slices);

movmax_val = ceil(1e3 / (hz));
moving_data = movmean(abs(chB), movmax_val);

new_data = reshape(moving_data, slice_size,[]);

% summerizing data -> making number repeating alot of times repeat once
all_slices = ones(1, number_of_slices);
for i=1:1:number_of_slices
    all_slices(i) = summerize_data(new_data(:,i));
end

% make the summerized data equal 0, 0.5 or 1
min_val = min(all_slices);
max_val = max(all_slices);
for i=1:1:number_of_slices
    all_slices(i) = discrete(all_slices(i), min_val, max_val);
end

% find first 0.5 bit:
first_half_range = find(all_slices == 0.5);
first_half = find_half(first_half_range, slice_per_bit);

save_indexes = ones(2, bit_amount - 1);

final_data = ones(1, bit_amount - 1);
for i=1:1:bit_amount - 1
    phase_shift = 0;
%     if i >= 22
%         phase_shift = 2;
%     end
%     if i >= 190
%         phase_shift = 0;
%     end
    first_index = mod(first_half + (i + phase_shift) * (slice_per_bit) - 1, number_of_slices) + 1;
    second_index =  mod(first_half + (i + 1 + phase_shift) * (slice_per_bit) - 1, number_of_slices) + 1;
%     if first_index == 0
%         first_index = 1;
%         second_index = second_index + 1;
%     end
    save_indexes(1, i) = first_index;
    save_indexes(2, i) = second_index;
    if first_index < second_index
    final_data(i) = mode(all_slices(first_index:second_index)');
    else
    array_to_check = cat(2, all_slices(first_index:number_of_slices), all_slices(1:second_index));
    final_data(i) = mode(array_to_check');
    end
end

% Channel A

axisHandleChA = subplot(3,1,1); 
hold on
plot(axisHandleChA, timeMs, chA, 'b');
% scatter(axisHandleChA, timeMs(save_indexes(1, :) * slice_size), zeros(1, length(save_indexes(1, :))), 'g*');
hold off

ylim(axisHandleChA, [-2500 2500]); % Adjust vertical axis for signal.
title(axisHandleChA, 'Channel A');
xlabel(axisHandleChA, 'Time (s)');
ylabel(axisHandleChA, 'Voltage (mV)');
grid(axisHandleChA);


% Channel B
axisHandleChB = subplot(3,1,2); 
hold on
% plot(axisHandleChB, timeMs, chB, 'r');
plot(axisHandleChB, timeMs, moving_data, 'r');
% scatter(axisHandleChB, timeMs(save_indexes(1, :) * slice_size), moving_data(save_indexes(1, :)), 'g*');
% ylim(axisHandleChB, [-1500 1500]); % Adjust vertical axis for signal.
hold off

title(axisHandleChB, 'Channel B');
xlabel(axisHandleChB, 'Time (s)');
ylabel(axisHandleChB, 'Voltage (mV)');
grid(axisHandleChB);

axisHandleBits = subplot(3,1,3);
hold on
plot(axisHandleBits, repmat(repelem(final_data, slice_size), 1, 1));
hold off

title(axisHandleBits, 'Decoding');
xlabel(axisHandleBits, 'Bit slices');
ylabel(axisHandleBits, 'Bit');
grid(axisHandleBits);


%% Turn off signal generator
% Sets the output to 0 V DC.

[status.setSigGenOff] = invoke(sigGenGroupObj, 'setSigGenOff');

%% Arbitrary waveform generator - output shots
% Output 2 cycles of an arbitrary waveform using a software trigger.
%
% Note that the signal generator will output the value coresponding to the
% first sample in the arbitrary waveform until the trigger event occurs.

% Increment      : 0 Hz
% Dwell Time     : 1 s
% Arb. Waveform  : y (defined above)
% Sweep Type     : 0 (ps2000aEnuminfo.enPS2000ASweepType.PS2000A_UP)
% Operation      : 0 (ps2000aEnuminfo.enPS2000AExtraOperations.PS2000A_ES_OFF)
% Index Mode     : 0 (ps2000aEnuminfo.enPS2000AIndexMode.PS2000A_SINGLE)
% Shots          : 2 
% Sweeps         : 0
% Trigger Type   : 0 (ps2000aEnuminfo.enPS2000ASigGenTrigType.PS2000A_SIGGEN_RISING)
% Trigger Source : 4 (ps2000aEnuminfo.enPS2000ASigGenTrigSource.PS2000A_SIGGEN_SOFT_TRIG)
% Ext. Threshold : 0 mV (Only supported by PicoScope 2206/7/8 models)

[status.setSigGenArbitrary] = invoke(sigGenGroupObj, 'setSigGenArbitrary', 1e-20, 1, y, 0, 0, 0, 1, 0, 0, 4, 0);

% Trigger the AWG

% State : 1 (a non-zero value will trigger the output)
[status.sigGenSoftwareControl] = invoke(sigGenGroupObj, 'ps2000aSigGenSoftwareControl', 1);

%% Turn off signal generator
% Sets the output to 0 V DC.

[status.setSigGenOff] = invoke(sigGenGroupObj, 'setSigGenOff');

%% Disconnect device
% Disconnect device object from hardware.

disconnect(ps2000aDeviceObj);
delete(ps2000aDeviceObj);
%%
% save data: 
distance = string(5);
writematrix(chA, "data\saved_data_2240_A_" + distance + ".csv");
writematrix(chB, "data\saved_data_2240_B_"+ distance + ".csv");

%% PLOTS: 

% hold on
% plot(final_data);
% plot(bits);
% hold off

%% 
axisHandleData = subplot(2,1,1);
hold on
image(axisHandleData ,flip(reshape(cat(2, [1] ,bits), 16, []), 1),'CDataMapping','scaled');
colormap("bone")
ylim(axisHandleData, [0.5 16.5]);
xlim(axisHandleData, [0.5 16.5]);
title(axisHandleData, 'Encoded Picture');
hold off

axisHandleFinalData = subplot(2,1,2);
for i=1:1:length(final_data)
    if final_data(i) == 0.5
        final_data(i) = 0;
    end
end
hold on
image(axisHandleFinalData ,flip(reshape(cat(2, [1] ,final_data), 16, []), 1),'CDataMapping','scaled');
ylim(axisHandleFinalData, [0.5 16.5]);
xlim(axisHandleFinalData, [0.5 16.5]);
hold off
title(axisHandleFinalData, 'Decoded Picture distance = ' + distance + " {cm}");

%%
function half = find_half(half_array, slice_per_bit)
    half = half_array(1);
    accuracy_factor = slice_per_bit / 2;
    for i=1:1:length(half_array) - accuracy_factor
        counter = 0;
        for j=1:1:accuracy_factor
            if half_array(i + j) - half_array(i + j - 1) == 1
                counter = counter + 1;
            end
        end
        if counter == accuracy_factor
            half = half_array(i);
            break
        end
    end
end


function d = discrete(value, min_val, max_val)
    value = (value - min_val) / (max_val - min_val);
    if value > 0.66
        d = 1;
    elseif value < 0.35
        d = 0;
    else
        d = 0.5;
    end
end

function x = summerize_data(data_array)
    data_array = abs(data_array);
    x = max(data_array);
end