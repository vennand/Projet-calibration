function [model, data] = GenerateModel_OneMarker(data)
% version 1

% File extracted from DoCi.s2mMod

%%%
% IMPORTANT
% Changed order of rotation at the arm from yxz to xyz, because of the way
% spatial_v2 interprets them.
%%%

% Informations générales
% root_actuated	1
% external_forces	0

import casadi.*

model.name = 'GravityOptimisation';

model.NB = 3;

model.jtype = {'Px', 'Py', 'Pz'};
model.parent = [0,1,2];

model.Xtree = {eye(6), eye(6), eye(6)};
model.I = {zeros(6,6),zeros(6,6),mcI(0.1,[0.0000000000 0.0000000000 0.0000000000],[0.0000000000 0.0000000000 0.0000000000;0.0000000000 0.0000000000 0.0000000000;0.0000000000 0.0000000000 0.0000000000])};

model.markers.name = {'WeightedLine1', 'WeightedLine2', 'WeightedLine3', 'WeightedLine4', 'CalibrationWand1', 'CalibrationWand2', 'CalibrationWand3', 'CalibrationWand4', 'CalibrationWand5', 'FreeFallMarker1'};

model.markers.parent = [3 3 3];
model.markers.coordinates = [0 0 0];

model.appearance.base = { 'line', [1.1 0 0; 0 0 0; 0 1.1 0; 0 0 0; 0 0 1.1]};

model.appearance.body{1} = {};
model.appearance.body{2} = {};
model.appearance.body{3} = {'sphere', model.markers.coordinates, 0.01};

model.nq = model.NB;
model.nx = model.nq+model.nq;

data.x0 = zeros(model.nx,1);

model.idx_q = 1:model.nq;
model.idx_v = model.nq+1:2*model.nq;

qmin_base = [-inf,-inf,-inf];
qmax_base = [ inf, inf, inf];

qdotmin_base = [-inf,-inf,-inf];
qdotmax_base = [ inf, inf, inf];

model.xmin = [qmin_base'; ... % q
              qdotmin_base']; % qdot
model.xmax = [qmax_base'; ... % q
              qdotmax_base']; % qdot

if isfield(data, 'gravity')
    model.gravity = data.gravity;
end
end
