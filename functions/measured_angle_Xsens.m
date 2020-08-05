function data = measured_angle_Xsens(data)

list_of_files = {'MT_0120050D-000-000_00B41A46.txt', 'MT_0120050D-001-000.txt', 'MT_0120050D-002-000.txt', 'MT_0120050D-003-000.txt', 'MT_0120050D-004-000.txt', 'MT_0120050D-005-000.txt', 'MT_0120050D-006-000.txt', 'MT_0120050D-007-000.txt', 'MT_0120050D-008-000.txt'};

switch data.angle_measured
    case 0
        filename = list_of_files{1};
    case 1
        data.angle_Xsens = 'NaN';
        data.angle_Xsens_corrected = 'NaN';
        data.gravity_Xsens = 'NaN';
        return
    case 2
        data.angle_Xsens = 'NaN';
        data.angle_Xsens_corrected = 'NaN';
        data.gravity_Xsens = 'NaN';
        return
    case 3
        filename = list_of_files{2};
    case 4
        filename = list_of_files{3};
    case 5
        filename = list_of_files{4};
    case 6
        filename = list_of_files{5};
    case 7
        filename = list_of_files{6};
    case 8
        filename = list_of_files{7};
    case 9
        filename = list_of_files{8};
    case 10
        filename = list_of_files{9};
end

fileID = fopen(['/home/andre/Projet calibration André/Projet calibration André/2020-07-24/' filename]);
Xsens_data = textscan(fileID, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','HeaderLines',13);
fclose(fileID);

% Acc_X  Acc_Y  Acc_Z
acceleration_data = cell2mat(Xsens_data(2:4));
average_acceleration = mean(acceleration_data);

if data.angle_measured == 0
    correction_rotation_matrix = vrrotvec2mat(vrrotvec(average_acceleration, data.gravity));
    save('../correction_rotation_matrix.mat','correction_rotation_matrix')
else
    load('../correction_rotation_matrix.mat','correction_rotation_matrix')
end
average_acceleration_corrected = correction_rotation_matrix*average_acceleration';

angle_deviation = angle_between_vectors(average_acceleration, data.gravity');
angle_deviation_corrected = angle_between_vectors(average_acceleration_corrected, data.gravity');

data.angle_Xsens = angle_deviation;
data.angle_Xsens_corrected = angle_deviation_corrected;
data.gravity_Xsens = norm(average_acceleration);
end