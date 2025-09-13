% Flow maps
Ac_1 = randn(2);
Ac_2 = randn(2);

% Flow set angles
flow_set_start_angle_1 = 0;
flow_set_end_angle_1   = pi;
flow_set_start_angle_2 = 0;
flow_set_end_angle_2   = pi;

flow_set_angles_1 = generateFlowBasedConeAngles2D(...
  flowMapMatrix=Ac_1, ...
  flowSetStartAngle=flow_set_start_angle_1, ...
  flowSetEndAngle=flow_set_end_angle_1)

flow_set_angles_2 = generateFlowBasedConeAngles2D(...
  flowMapMatrix=Ac_2, ...
  flowSetStartAngle=flow_set_start_angle_2, ...
  flowSetEndAngle=flow_set_end_angle_2)

% Discrete jump maps.
Ad_1_mapsto_1 = rand(2); 
Ad_1_mapsto_2 = eye(2); 
Ad_2_mapsto_1 = rand(2); 
Ad_2_mapsto_2 = eye(2); 

jump_set_angles_1_mapsto_1 = [pi-0.1, pi];
jump_set_angles_1_mapsto_2 = [pi-0.1, pi];
jump_set_angles_2_mapsto_1 = [pi-0.1, pi];
jump_set_angles_2_mapsto_2 = [pi-0.1, pi];

jump_set_image_angles_1_mapsto_1 = mapAnglesByLinearMap(Ad_1_mapsto_1, jump_set_angles_1_mapsto_1)
jump_set_image_angles_1_mapsto_2 = mapAnglesByLinearMap(Ad_1_mapsto_2, jump_set_angles_1_mapsto_2)
jump_set_image_angles_2_mapsto_1 = mapAnglesByLinearMap(Ad_2_mapsto_1, jump_set_angles_2_mapsto_1)
jump_set_image_angles_2_mapsto_2 = mapAnglesByLinearMap(Ad_2_mapsto_2, jump_set_angles_2_mapsto_2)

mode_1_angles = [ ... 
  flow_set_angles_1, ... 
  jump_set_angles_1_mapsto_1, jump_set_angles_1_mapsto_2, ... 
  jump_set_image_angles_1_mapsto_1, jump_set_image_angles_2_mapsto_1 ... 
]
mode_2_angles = [ ... 
  flow_set_angles_2, ... 
  jump_set_angles_2_mapsto_1, jump_set_angles_2_mapsto_2, ... 
  jump_set_image_angles_1_mapsto_2, jump_set_image_angles_2_mapsto_2 ... 
] 

all_angles = [mode_1_angles, mode_2_angles];

conical_partition_1 = ConicalPartition.fromAngles2D(mode_1_angles)
conical_partition_2 = ConicalPartition.fromAngles2D(mode_2_angles)
conical_partition = ConicalPartition.fromAngles2D(all_angles);

figure(1);
clf();
xlim(1.0*[-1, 1]);
ylim(1.0*[-1, 1]);
axis square;
hold on;
conical_partition_1.plot()


figure(2);
clf();
xlim(1.0*[-1, 1]);
ylim(1.0*[-1, 1]);
axis square;
hold on;
conical_partition_2.plot()

figure(3);
clf();
xlim(1.0*[-1, 1]);
ylim(1.0*[-1, 1]);
axis square;
hold on;
conical_partition.plot()


conic_abstraction_1 = ConicAbstraction(conical_partition_1, flow_set_cone_ndxs, jump_set_cone_ndxs, jump_set_image_cone_ndxs, A_c, A_d);

warning("Returning early before the end of run_conic_abstraction_with_modes.");
return % TODO: Remove

conic_abstraction_1_mapsto_1 = ConicAbstraction.fromAngles2D(...
  flowMapMatrix=-eye(2), ...
  jumpMapMatrix=Ad_1_mapsto_1, ...
  jumpSetAngles=jump_set_angles_1_mapsto_1, ...
  flowSetAngles=flow_set_angles_1)
conic_abstraction_1_mapsto_2 = ConicAbstraction.fromAngles2D(...
  flowMapMatrix=-eye(2), ...
  jumpMapMatrix=Ad_1_mapsto_2, ...
  jumpSetAngles=jump_set_angles_1_mapsto_2, ...
  flowSetAngles=flow_set_angles_1)
conic_abstraction_2_mapsto_1 = ConicAbstraction.fromAngles2D(...
  flowMapMatrix=-eye(2), ...
  jumpMapMatrix=Ad_2_mapsto_1, ...
  jumpSetAngles=jump_set_angles_2_mapsto_1, ...
  flowSetAngles=flow_set_angles_2)
conic_abstraction_2_mapsto_2 = ConicAbstraction.fromAngles2D(...
  flowMapMatrix=-eye(2), ...
  jumpMapMatrix=Ad_2_mapsto_2, ...
  jumpSetAngles=jump_set_angles_2_mapsto_2, ...
  flowSetAngles=flow_set_angles_2)

% function 