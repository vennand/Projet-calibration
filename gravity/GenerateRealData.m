function [model, data] = GenerateRealData(model,data)
% frames = data.frames;
labels = data.labels;

real_data = ezc3dRead(data.dataFile);
frequency = real_data.header.points.frameRate;

one_marker = real_data.data.points(1,labels(end),:)/1000; %meter
first_frame = find(~isnan(one_marker), 1, 'first');
last_frame = find(~isnan(one_marker), 1, 'last');

frames = (first_frame+5):last_frame;
data.frames = frames;

data.realNint = length(data.frames);
data = adjust_number_of_interval(data);

data.Duration = length(frames)/frequency;

markers = real_data.data.points(:,labels,frames)/1000; %meter
labels_name = real_data.parameters.POINT.LABELS.DATA(labels);

% Nint = data.Nint;
% [num_dimension, num_label, num_frames] = size(markers);
[~, order] = ismember(model.markers.name, labels_name);

%Reorder the labels according to the model and to number of selected nodes
markers_reformat = markers(:, order, 1:data.step:end);

%Reposition the referential to the model
% ref_rotation = refential_matrix()';
% for k=1:size(markers_reformat,3)
%     markers_reformat(:,:,k) = ref_rotation * markers_reformat(:,:,k);
% end

% Storing data
data.markers = markers_reformat;
end
