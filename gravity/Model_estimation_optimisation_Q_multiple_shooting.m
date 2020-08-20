% Script to optimize a trajectory with 42 DoF, 1sec time frame
% models with trapezoidal collocation
clear, clc, close all
tic
run('../startup.m')
import casadi.*

data.nDoF = 3;

data.Nint = 30;% number of control nodes
data.odeMethod = 'rk4';
data.NLPMethod = 'MultipleShooting';

data.gravity = [0; 0; -9.81];
% data.gravity = [0; 0; -9.80639]; % According to WolframAlpha
data.gravityRotationBound = pi/16;
data.nCardinalCoor = 3;

output_filename = 'Solutions/angle_and_gravity.xls';
header = {'Trial', 'Angle smartphone', 'Angle fil à plomb', 'Angle Xsens', 'Angle Xsens corrigé', 'Angle gravité optimisée contrainte stricte', 'Angle gravité optimisée sans contrainte', 'Norme gravité Xsens', 'Norme gravité contrainte stricte', 'Norme gravité sans contrainte'};
writecell(header, output_filename)

% data.angle_measured = 0;
% data.trial_type = 'Para';
% data.trial_number = 2;

% angles_measured = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
angles_measured = 10;
trial_types = {'Drop', 'Para'};
trial_numbers = [1, 2, 3, 4];
trial_to_ignore = {'Angle00_Para01', 'Angle00_Para04', 'Angle01_Drop01', 'Angle02_Para02', 'Angle02_Para03', ...
                   'Angle03_Para04', 'Angle04_Para03', 'Angle04_Para04', 'Angle06_Para04', 'Angle09_Para04'};

constrained_bool = [true, false];

for constrained = constrained_bool
for angle_measured = angles_measured
    data.angle_measured = angle_measured;
for trial_type_cell = trial_types
    trial_type = trial_type_cell{1};
    data.trial_type = trial_type;
for trial_number = trial_numbers
    data.trial_number = trial_number;
    if strcmp(data.trial_type, 'Drop') && data.trial_number == 4
        continue
    end

    if strcmp(data.trial_type, 'Para')
        increment_trial_type = 3;
    else
        increment_trial_type = 0;
    end
    index = 2 + 7*(data.angle_measured) + increment_trial_type + (data.trial_number-1);
    
    file = ['Angle' num2str(data.angle_measured,'%02d') '_' data.trial_type num2str(data.trial_number,'%02d')];
    
    if any(strcmp(trial_to_ignore, file))
        output_file = readcell(output_filename);
        output_file{index,1} = file;
        output_file(index,2:end) = {'NaN'};
        writecell(output_file, output_filename)
        continue
    end

data.dataFile = ['../Projet calibration André/2020-07-24/Calibration Mocap/New Patient/New Session/' file '.c3d'];
data.kalmanDataFile_q = ['../EKF/Angle' num2str(data.angle_measured,'%02d') '_' data.trial_type num2str(data.trial_number,'%02d') '_Q.mat'];
data.kalmanDataFile_v = ['../EKF/Angle' num2str(data.angle_measured,'%02d') '_' data.trial_type num2str(data.trial_number,'%02d') '_Qd.mat'];
data.kalmanDataFile_a = ['../EKF/Angle' num2str(data.angle_measured,'%02d') '_' data.trial_type num2str(data.trial_number,'%02d') '_Qdd.mat'];

% Spécific à Do_822_contact_2.c3d
% Le drop est entre les frames 81 et 142, mais Kalman n'est pas bon à 81
% data.frames = 82:142;
data.labels = 1:10;

data.weightU = 1e-7;
data.weightX = 1;
data.weightQV = [1; 0.01];
data.weightPoints = 1;

disp('Generating Model')
[model, data] = GenerateModel_OneMarker(data);
disp('Loading Real Data')
[model, data] = GenerateRealData(model,data);
disp('Loading Kalman Filter')
[model, data] = GenerateKalmanFilter(model,data);
disp('Calculating Estimation')
if constrained
    [prob, lbw, ubw, lbg, ubg, objFunc, conFunc, objGrad, conGrad] = GenerateEstimation_Q_multiple_shooting(model, data);
else
    [prob, lbw, ubw, lbg, ubg, objFunc, conFunc, objGrad, conGrad] = GenerateEstimation_Q_multiple_shooting_unconstrained(model, data);
end

% [lbw, ubw] = GenerateInitialConstraints(model, data, lbw, ubw);
% [lbw, ubw] = GenerateFinalConstraints(model, data, lbw, ubw);

options = struct;
options.ipopt.max_iter = 3000;
options.ipopt.print_level = 5;
options.ipopt.linear_solver = 'ma57';

options.ipopt.tol = 1e-6; % default: 1e-08
% options.ipopt.acceptable_tol = 1e-4; % default: 1e-06
options.ipopt.constr_viol_tol = 0.001; % default: 0.0001
% options.ipopt.acceptable_constr_viol_tol = 0.1; % default: 0.01

disp('Generating Solver')
% solver = nlpsol('solver', 'snopt', prob, options); % FAIRE MARCHER ÇA
solver = nlpsol('solver', 'ipopt', prob, options);

w0=[];
for k=1:data.Nint+1
    w0 = [w0; data.kalman_q(:,k); data.kalman_v(:,k)];
end

N_G = data.nCardinalCoor;
w0 = [w0; data.gravity];

sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, 'lbg', lbg, 'ubg', ubg);

q_opt = nan(model.nq,data.Nint+1);
v_opt = nan(model.nq,data.Nint+1);
w_opt = full(sol.x);

for i=1:model.nq
    q_opt(i,:) = w_opt(i:model.nx:end - N_G)';
    v_opt(i,:) = w_opt(i+model.nq:model.nx:end - N_G)';
end
G_opt = w_opt(end - N_G + 1:end);

data.G_opt = G_opt;

data.q_opt = q_opt;
data.v_opt = v_opt;

model.gravity = data.G_opt;

% disp('Calculating Simulation')
% [model, data] = GenerateSimulation(model, data);
disp('Calculating Momentum')
data = CalculateMomentum(model, data);

stats = solver.stats;

angle_deviation = angle_between_vectors(data.G_opt, data.gravity);
data.angle_deviation = angle_deviation;

data = measure_angle_weighted_line(data);
data = measured_angle_Xsens(data);

save(['Solutions/Angle' num2str(data.angle_measured,'%02d') '_'  data.trial_type num2str(data.trial_number,'%02d') '_Constrainted_' num2str(constrained) '.mat'],'model','data','stats')
% GeneratePlots(model, data);
% AnimatePlot(model, data, 'sol', 'kalman');

% disp('Angle optimisé de déviation de la gravité')
% disp([num2str(data.angle_deviation, '%1.5f') ' degrés'])
% 
% disp('Angle du fil à pèche pesé par rapport à la calibration')
% disp([num2str(data.angle_weighted_line, '%1.5f') ' degrés'])
% 
% disp('Angle de la baguette mesuré par le Xsens')
% disp([num2str(data.angle_Xsens, '%1.5f') ' degrés'])
% 
% disp('Angle de la baguette mesuré par le Xsens, corrigé')
% disp([num2str(data.angle_Xsens_corrected, '%1.5f') ' degrés'])
% 
% disp('Angle de la baguette mesuré par le smartphone')
% disp([num2str(data.angle_measured) ' degrés'])
% 
% disp('Norme de la gravité optimisée')
% disp(norm(data.G_opt))
% 
% disp('Norme de la gravité Xsens')
% disp(norm(data.gravity_Xsens))

output_file = readcell(output_filename);

if any(strcmp(output_file(:,1), file))
    if constrained
        output_values1 = {data.angle_measured, data.angle_weighted_line, data.angle_Xsens, data.angle_Xsens_corrected, data.angle_deviation};
        output_values2 = {data.gravity_Xsens, norm(data.G_opt)};

        output_file(index,2:6) = output_values1;
        output_file(index,8:9) = output_values2;
    else
        output_values1 = {data.angle_deviation};
        output_values2 = {norm(data.G_opt)};

        output_file(index,7) = output_values1;
        output_file(index,10) = output_values2;
    end
else
    output_file{index,1} = file;
    if constrained
        output_values1 = {data.angle_measured, data.angle_weighted_line, data.angle_Xsens, data.angle_Xsens_corrected, data.angle_deviation};
        output_values2 = {data.gravity_Xsens, norm(data.G_opt)};

        output_file(index,2:6) = output_values1;
        output_file(index,7) = {'NaN'};
        output_file(index,8:9) = output_values2;
        output_file(index,10) = {'NaN'};
    else
        output_values1 = {data.angle_deviation};
        output_values2 = {norm(data.G_opt)};

        output_file(index,2:6) = {'NaN'};
        output_file(index,7) = output_values1;
        output_file(index,8:9) = {'NaN'};
        output_file(index,10) = output_values2;
    end
end

writecell(output_file, output_filename)
disp(file)

end
end
end
end
toc
% showmotion(model, 0:data.Duration/data.Nint:data.Duration, data.q_opt(:,:))
