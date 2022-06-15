% hold on
% plot(reshape(input, [], 1));
% plot(reshape(output, [], 1));
% hold off
%%
imdata = imread('ok.bmp');
input = double(imdata);
bit_amount = 256;
input = reshape(input, [], bit_amount);
inputLen = 1024;
output = ones(1, inputLen);
output = reshape(output, [] ,bit_amount);

clc;
close all;

PS2000aConfig;

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


for i=1:1:(inputLen / bit_amount)
    output(i, :) = recieveBits(input(i, :), ps2000aDeviceObj,ps2000aEnuminfo);
end

disconnect(ps2000aDeviceObj);
delete(ps2000aDeviceObj);

%% 
edge_size = sqrt(inputLen);

axisHandleData = subplot(2,1,1);
hold on
image(axisHandleData ,flip(reshape(input, edge_size, []), 1),'CDataMapping','scaled');
colormap("bone")
ylim(axisHandleData, [0.5 (edge_size +0.5 )]);
xlim(axisHandleData, [0.5 (edge_size + 0.5)]);
title(axisHandleData, 'Encoded Picture');
hold off

axisHandleFinalData = subplot(2,1,2);
for i=1:1:length(output)
    if output(i) == 0.5
        output(i) = 0;
    end
end

hold on
image(axisHandleFinalData ,flip(reshape(output, edge_size, []), 1),'CDataMapping','scaled');
ylim(axisHandleFinalData, [0.5 (edge_size + 0.5)]);
xlim(axisHandleFinalData, [0.5 (edge_size + 0.5)]);
hold off
title(axisHandleFinalData, 'Decoded Picture = ' + string(hz));


%%
function outputBits = recieveBits(inputBits, ps2000aDeviceObj, ...
    ps2000aEnuminfo)
sigGenGroupObj = get(ps2000aDeviceObj, 'Signalgenerator');
sigGenGroupObj = sigGenGroupObj(1);

[status.setSigGenOff] = invoke(sigGenGroupObj, 'setSigGenOff');

% Configure property value(s).
plaster_freq = 1.024;
hz = 1;
set(ps2000aDeviceObj.Signalgenerator(1), 'startFrequency', 1.0 / plaster_freq * hz);
set(ps2000aDeviceObj.Signalgenerator(1), 'stopFrequency', 1.0 / plaster_freq * hz);
set(ps2000aDeviceObj.Signalgenerator(1), 'offsetVoltage', 0.0);
set(ps2000aDeviceObj.Signalgenerator(1), 'peakToPeakVoltage', 2000.0);

awgBufferSize = get(sigGenGroupObj, 'awgBufferSize');

% pause(0.5);
bit_amount = 256;

x = linspace(0, 2 * pi, awgBufferSize / bit_amount);
zero_bit = zeros(1, awgBufferSize / bit_amount);
frequency = 341.4 * 1e3;
one_bit = sin(x * frequency) * 2;

bits = inputBits(2:bit_amount);
% bits = randi([0 1],1, bit_amount - 1);
y = [];
y = cat(2, y, one_bit / 3);
for bit=bits
    if bit == 0
        y = cat(2, y, zero_bit);
    end
    if bit == 1
        y = cat(2, y, one_bit);
    end
end


[status.setSigGenArbitrarySimple] = invoke(sigGenGroupObj, 'setSigGenArbitrarySimple', y);

blockGroupObj = get(ps2000aDeviceObj, 'Block');
blockGroupObj = blockGroupObj(1);
set(ps2000aDeviceObj, 'numPostTriggerSamples', 8192 * 125 / hz);

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
    save_indexes(1, i) = first_index;
    save_indexes(2, i) = second_index;
    if first_index < second_index
    final_data(i) = mode(all_slices(first_index:second_index)');
    else
    array_to_check = cat(2, all_slices(first_index:number_of_slices), all_slices(1:second_index));
    final_data(i) = mode(array_to_check');
    end
end
for i=1:1:length(final_data)
    if final_data(i) == 0.5
        final_data(i) = 0;
    end
end
outputBits = cat(2, [1], final_data);
% Channel A

% axisHandleChA = subplot(3,1,1); 
% hold on
% plot(axisHandleChA, timeMs, chA, 'b');
% % scatter(axisHandleChA, timeMs(save_indexes(1, :) * slice_size), zeros(1, length(save_indexes(1, :))), 'g*');
% hold off
% 
% ylim(axisHandleChA, [-2500 2500]); % Adjust vertical axis for signal.
% title(axisHandleChA, 'Channel A');
% xlabel(axisHandleChA, 'Time (s)');
% ylabel(axisHandleChA, 'Voltage (mV)');
% grid(axisHandleChA);
% 
% 
% % Channel B
% axisHandleChB = subplot(3,1,2); 
% hold on
% % plot(axisHandleChB, timeMs, chB, 'r');
% plot(axisHandleChB, timeMs, moving_data, 'r');
% % scatter(axisHandleChB, timeMs(save_indexes(1, :) * slice_size), moving_data(save_indexes(1, :)), 'g*');
% % ylim(axisHandleChB, [-1500 1500]); % Adjust vertical axis for signal.
% hold off
% 
% title(axisHandleChB, 'Channel B');
% xlabel(axisHandleChB, 'Time (s)');
% ylabel(axisHandleChB, 'Voltage (mV)');
% grid(axisHandleChB);
% 
% axisHandleBits = subplot(3,1,3);
% hold on
% plot(axisHandleBits, repmat(repelem(final_data, slice_size), 1, 1));
% hold off
% 
% title(axisHandleBits, 'Decoding');
% xlabel(axisHandleBits, 'Bit slices');
% ylabel(axisHandleBits, 'Bit');
% grid(axisHandleBits);


[status.setSigGenOff] = invoke(sigGenGroupObj, 'setSigGenOff');

[status.setSigGenArbitrary] = invoke(sigGenGroupObj, 'setSigGenArbitrary', 1e-20, 1, y, 0, 0, 0, 1, 0, 0, 4, 0);

[status.sigGenSoftwareControl] = invoke(sigGenGroupObj, 'ps2000aSigGenSoftwareControl', 1);
end

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