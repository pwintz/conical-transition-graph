
clear;

% Flow maps
% mat_1_eig_vecs = [
%   [1; -1], []
% ]
% mat_2_eig_vals = [
%   1.0, 
%   -2
% ]
flow_map_matrix_1 = [
   2, 2;
  -3, 1;
];
flow_map_matrix_2 = [
  -1, 1;
  -4, -2;
];
f1 = @(x) flow_map_matrix_1 * x;
f2 = @(x) flow_map_matrix_2 * x;

[eigenvectors_1, eigenvalues_1] = eig(flow_map_matrix_1);
[eigenvectors_2, eigenvalues_2] = eig(flow_map_matrix_2);
eigvec_ndx_in_right_half_plane = eigenvectors_2(:, diag(eigenvalues_1) > 0);
% pwintz.assertions.assertIsScalar(eigvec_ndx_in_right_half_plane);
% angle_to_avoid = pwintz.math.atan2(eigenvectors_2(:, eigvec_ndx_in_right_half_plane))

% ╭─────────────────────────────────────────╮
% │             Flow set angles             │
% ╰─────────────────────────────────────────╯
flow_set_start_angle_1 = pi/2; flow_set_end_angle_1   = 3*pi / 2;
flow_set_start_angle_2 = 0;    flow_set_end_angle_2   = 2*pi;

flow_set_angles_1 = generateFlowBasedConeAngles2D(...
  flowMapMatrix=flow_map_matrix_1, ...
  flowSetStartAngle=flow_set_start_angle_1, ...
  flowSetEndAngle=flow_set_end_angle_1, ...
  maxStateConeAngle = 2*pi/10, ...
  maxDerivativeConeAngle = 2*pi/30 ...
);

flow_set_angles_2 = generateFlowBasedConeAngles2D(...
  flowMapMatrix=flow_map_matrix_2, ...
  flowSetStartAngle=flow_set_start_angle_2, ...
  flowSetEndAngle=flow_set_end_angle_2,...
  maxStateConeAngle = 2*pi/10, ...
  maxDerivativeConeAngle = 2*pi/30 ...
);

% ⋘──────── Jump image angles ────────⋙
jump_set_1_mapsto_1_start_angle  = mod(pwintz.math.atan2([ -1;  0]), 2*pi);
jump_set_1_mapsto_1_end_angle    = mod(pwintz.math.atan2([ -4; -1]), 2*pi);

jump_set_1_mapsto_2_start_angle  = mod(pwintz.math.atan2([  0;  1]), 2*pi);
jump_set_1_mapsto_2_end_angle    = mod(pwintz.math.atan2([ -1;  2]), 2*pi);

jump_set_2_mapsto_1_start_angle  = mod(pwintz.math.atan2([ -1;  0]), 2*pi);
jump_set_2_mapsto_1_end_angle    = mod(pwintz.math.atan2([ -4; -1]), 2*pi);

jump_set_2_mapsto_2_start_angle  = mod(pwintz.math.atan2([  1;  0]), 2*pi);
jump_set_2_mapsto_2_end_angle    = mod(pwintz.math.atan2([  4;  1]), 2*pi);

% pwintz.math.atan2([-1; 0]);  pi;         
% pwintz.math.atan2([-4; -1]); pi + pi/12;
% pwintz.math.atan2([0; 1]);   pi/2;      
% pwintz.math.atan2([-1; 2]);  2*pi/3;
% pwintz.math.atan2([-4; 1]);  pi - pi/12; 
% pwintz.math.atan2([-1; 0]);  pi;
% pwintz.math.atan2([1; 0]);   0;          
% pwintz.math.atan2([4; 1]);   pi/12;

% Check angles
pwintz.assertions.assertAllLessThan(jump_set_1_mapsto_1_start_angle, jump_set_1_mapsto_1_end_angle);
pwintz.assertions.assertAllLessThan(jump_set_1_mapsto_2_start_angle, jump_set_1_mapsto_2_end_angle);
pwintz.assertions.assertAllLessThan(jump_set_2_mapsto_1_start_angle, jump_set_2_mapsto_1_end_angle);
pwintz.assertions.assertAllLessThan(jump_set_2_mapsto_2_start_angle, jump_set_2_mapsto_2_end_angle);


gammas = 10.^(-3:0.5:1)

gamma_results = dictionary();

for gamma = gammas
  
  % ╭───────────────────────────────────────╮
  % │             Jump Matrices             │
  % ╰───────────────────────────────────────╯
  jump_map_matrix_1_mapsto_1 = [2, 1/2; -2, 2]; % [-1.6, 3; -0.2, -0.1]; 
  jump_map_matrix_1_mapsto_2 = gamma*[1, 1; 0, 1]; 
  jump_map_matrix_2_mapsto_1 = [1, 3; 4, 2]; 
  jump_map_matrix_2_mapsto_2 = pwintz.linear_algebra.rotation2(-pi/2); 

  % ⋘──────── Jump set image angles ────────⋙
  jump_set_image_angles_1_mapsto_1 = mapAnglesByLinearMap(jump_map_matrix_1_mapsto_1, [jump_set_1_mapsto_1_start_angle, jump_set_1_mapsto_1_end_angle]);
  jump_set_image_angles_1_mapsto_2 = mapAnglesByLinearMap(jump_map_matrix_1_mapsto_2, [jump_set_1_mapsto_2_start_angle, jump_set_1_mapsto_2_end_angle]);
  jump_set_image_angles_2_mapsto_1 = mapAnglesByLinearMap(jump_map_matrix_2_mapsto_1, [jump_set_2_mapsto_1_start_angle, jump_set_2_mapsto_1_end_angle]);
  jump_set_image_angles_2_mapsto_2 = mapAnglesByLinearMap(jump_map_matrix_2_mapsto_2, [jump_set_2_mapsto_2_start_angle, jump_set_2_mapsto_2_end_angle]);

  mode_1_angles = [ ... 
    ... % Add angles for flows.
    flow_set_angles_1, ... 
    ... % Add angles for jumps from this mode.
    jump_set_1_mapsto_1_start_angle, ...
    jump_set_1_mapsto_1_end_angle, ...
    jump_set_1_mapsto_2_start_angle, ...
    jump_set_1_mapsto_2_end_angle, ...
    ... % Add angles for jumps into this mode.
    jump_set_image_angles_1_mapsto_1, ... 
    jump_set_image_angles_2_mapsto_1 ... 
  ];

  mode_2_angles = [ ... 
    ... % Add angles for flows.
    flow_set_angles_2, ... 
    ... % Add angles for jumps from this mode.
    jump_set_2_mapsto_1_start_angle, ...
    jump_set_2_mapsto_1_end_angle, ...
    jump_set_2_mapsto_2_start_angle, ...
    jump_set_2_mapsto_2_end_angle, ...
    ... % Add angles for jumps into this mode.  
    jump_set_image_angles_1_mapsto_2, ... 
    jump_set_image_angles_2_mapsto_2 ... 
  ];

  all_angles = [mode_1_angles, mode_2_angles];

  conical_partition_1 = ConicalPartition.fromAngles2D(mode_1_angles);
  conical_partition_2 = ConicalPartition.fromAngles2D(mode_2_angles);

  % flow_set_cone_ndxs_1 = conical_partition_1.cone_indices;
  % flow_set_cone_ndxs_2 = find(flow_set_angles_2 <= pi);


  % ⋘──────── Flow set indices ────────⋙
  flow_set_1_cone_ndxs = conical_partition_1.getConesIntersectingArc(flow_set_start_angle_1, flow_set_end_angle_1);
  flow_set_2_cone_ndxs = conical_partition_2.getConesIntersectingArc(flow_set_start_angle_2, flow_set_end_angle_2);

  % ⋘──────── Jump set indices (starts of jumps) ────────⋙
  jump_set_1_mapsto_1_cone_ndxs = conical_partition_1.getConesIntersectingArc(jump_set_1_mapsto_1_start_angle , jump_set_1_mapsto_1_end_angle);
  jump_set_1_mapsto_2_cone_ndxs = conical_partition_1.getConesIntersectingArc(jump_set_1_mapsto_2_start_angle , jump_set_1_mapsto_2_end_angle);
  jump_set_2_mapsto_1_cone_ndxs = conical_partition_2.getConesIntersectingArc(jump_set_2_mapsto_1_start_angle , jump_set_2_mapsto_1_end_angle);
  jump_set_2_mapsto_2_cone_ndxs = conical_partition_2.getConesIntersectingArc(jump_set_2_mapsto_2_start_angle , jump_set_2_mapsto_2_end_angle);

  pwintz.assertions.assertAllAreMembers(jump_set_1_mapsto_1_cone_ndxs, conical_partition_1.cone_indices);
  pwintz.assertions.assertAllAreMembers(jump_set_1_mapsto_2_cone_ndxs, conical_partition_1.cone_indices);
  pwintz.assertions.assertAllAreMembers(jump_set_2_mapsto_1_cone_ndxs, conical_partition_2.cone_indices);
  pwintz.assertions.assertAllAreMembers(jump_set_2_mapsto_2_cone_ndxs, conical_partition_2.cone_indices);

  pwintz.assertions.assertNonempty(jump_set_1_mapsto_1_cone_ndxs);
  pwintz.assertions.assertNonempty(jump_set_1_mapsto_2_cone_ndxs);
  pwintz.assertions.assertNonempty(jump_set_2_mapsto_1_cone_ndxs);
  pwintz.assertions.assertNonempty(jump_set_2_mapsto_2_cone_ndxs);

  % % Jump set image indices
  % jump_set_1_mapsto_1_image_cone_ndxs = conical_partition_1.getConesIntersectingArc(jump_set_1_mapsto_1_image_start_angle , jump_set_1_mapsto_1_image_end_angle);
  % jump_set_1_mapsto_2_image_cone_ndxs = conical_partition_1.getConesIntersectingArc(jump_set_1_mapsto_2_image_start_angle , jump_set_1_mapsto_2_image_end_angle);
  % jump_set_2_mapsto_1_image_cone_ndxs = conical_partition_1.getConesIntersectingArc(jump_set_2_mapsto_1_image_start_angle , jump_set_2_mapsto_1_image_end_angle);
  % jump_set_2_mapsto_2_image_cone_ndxs = conical_partition_1.getConesIntersectingArc(jump_set_2_mapsto_2_image_start_angle , jump_set_2_mapsto_2_image_end_angle);

  flow_specifications = [
    FlowSpecification(flow_set_1_cone_ndxs, flow_map_matrix_1);
    FlowSpecification(flow_set_2_cone_ndxs, flow_map_matrix_2);
  ];

  jump_specifications = [
    % First row defines jumps from mode 1.
    JumpSpecification(jump_set_1_mapsto_1_cone_ndxs, jump_map_matrix_1_mapsto_1), JumpSpecification(jump_set_1_mapsto_2_cone_ndxs, jump_map_matrix_1_mapsto_2);
    % Second row defines jumps from mode 2.
    JumpSpecification(jump_set_2_mapsto_1_cone_ndxs, jump_map_matrix_2_mapsto_1), JumpSpecification(jump_set_2_mapsto_2_cone_ndxs, jump_map_matrix_2_mapsto_2);
  ];

  conical_partitions = [
    conical_partition_1;
    conical_partition_2;
  ];

  jump_graph = JumpTransitionGainDigraph(conical_partitions, jump_specifications);
  flow_graph = FlowTransitionGainDigraph(conical_partition_1, flow_set_1_cone_ndxs, flow_map_matrix_1, conical_partition_2, flow_set_2_cone_ndxs, flow_map_matrix_2);


  jump_set_indices_in_mode_1 = [jump_set_1_mapsto_1_cone_ndxs, jump_set_1_mapsto_2_cone_ndxs];
  jump_set_indices_in_mode_2 = [jump_set_2_mapsto_1_cone_ndxs, jump_set_2_mapsto_2_cone_ndxs];
  jump_set_image_indices_in_mode_1 = [jump_graph.jump_set_image_cone_ndxs{:, 1}];
  jump_set_image_indices_in_mode_2 = [jump_graph.jump_set_image_cone_ndxs{:, 2}];

  flow_graph_between_jump_cones = flow_graph.contractEdgesToBetweenCones(...
    ... % v---Start indices ---v          v ---End indices ---v
    jump_set_image_indices_in_mode_1, jump_set_indices_in_mode_1, ...
    jump_set_image_indices_in_mode_2, jump_set_indices_in_mode_2 ...
  );

  
  %%
  ctg = TransitionGainDigraph.union(flow_graph_between_jump_cones, jump_graph);
  [cycle_min_gains, cycle_max_gains, cycles_nodes, cycles_edges] = ctg.getCycleGains();
  %%
  % pwintz.strings.format("There are %d cycles.", numel(cycle_min_gains));
  % pwintz.strings.format("minimum gain: %s", min(cycle_min_gains));
  % pwintz.strings.format("maximum gain: %s", max(cycle_max_gains));

  gamma_results{gamma} = struct("min_gain", min(cycle_min_gains), ...
                                "max_gain", max(cycle_max_gains), ...
                                "n_cycles", numel(cycle_min_gains));

  % gamma_results{gamma}            

  % pwintz.plots.namedFigure("Gains");
  % histogram(cycle_min_gains) 
  % end

end

return

%%
% close all; % Close figures.
pwintz.plots.namedFigure("Conical Transition Graph");
clf();
xlim(1.9*[-1, 1]);
ylim(0.8*[-1, 1]);
axis equal;
hold on;
axis off ;

% Turn off clipping
ax = gca;
ax.Clipping = "off";

% conical_partition_1.plot(originOff);
mode_1_plot_center = [-1.2; 0];
mode_2_plot_center = -mode_1_plot_center;
conical_partition_1.plot(offset=mode_1_plot_center);
conical_partition_2.plot(offset=mode_2_plot_center);

conical_partition_1.plotCone(flow_set_1_cone_ndxs, FaceColor="blue", FaceAlpha=0.2, offset = mode_1_plot_center);
conical_partition_2.plotCone(flow_set_2_cone_ndxs, FaceColor="blue", FaceAlpha=0.2, offset = mode_2_plot_center);
conical_partition_1.plotCone(jump_set_1_mapsto_1_cone_ndxs, FaceColor="red", FaceAlpha=0.2, offset = mode_1_plot_center);
conical_partition_1.plotCone(jump_set_1_mapsto_2_cone_ndxs, FaceColor="red", FaceAlpha=0.2, offset = mode_1_plot_center);
conical_partition_2.plotCone(jump_set_2_mapsto_1_cone_ndxs, FaceColor="red", FaceAlpha=0.2, offset = mode_2_plot_center);
conical_partition_2.plotCone(jump_set_2_mapsto_2_cone_ndxs, FaceColor="red", FaceAlpha=0.2, offset = mode_2_plot_center);


% jump_graph.jump_set_image_cone_ndxs{1, 1}
% jump_graph.jump_set_image_cone_ndxs{:, 2}

conical_partition_1.plotCone(jump_graph.jump_set_image_cone_ndxs{1, 1}, FaceColor="yellow", FaceAlpha=0.2, offset = mode_1_plot_center);
conical_partition_2.plotCone(jump_graph.jump_set_image_cone_ndxs{1, 2}, FaceColor="yellow", FaceAlpha=0.2, offset = mode_2_plot_center);
conical_partition_1.plotCone(jump_graph.jump_set_image_cone_ndxs{2, 1}, FaceColor="yellow", FaceAlpha=0.2, offset = mode_1_plot_center);
conical_partition_2.plotCone(jump_graph.jump_set_image_cone_ndxs{2, 2}, FaceColor="yellow", FaceAlpha=0.2, offset = mode_2_plot_center);

% flow_graph.plot();

pwintz.plots.plotLinearVectorField(flow_map_matrix_1, offset=mode_1_plot_center, xRange=1.1*[-1, 1], yRange=1.2*[-1, 1], Color=0.5*[1 1 1]);
pwintz.plots.plotLinearVectorField(flow_map_matrix_2, offset=mode_2_plot_center, xRange=1.1*[-1, 1], yRange=1.2*[-1, 1], Color=0.5*[1 1 1]);

jump_graph.plot("EdgeColor", "red");
flow_graph_between_jump_cones.plot("EdgeColor", "blue");
% ctg.plot();


pwintz.plots.saveExampleFigure('C:\Users\pwintz\code\Zeno\latex\images', 'ctg_two_modes', width=500, height=400);

%%

% pwintz.plots.namedFigure("Contracted Flow Graphs");
% clf();
% xlim(2.5*[-1, 1]);
% ylim(1.0*[-1, 1]);
% axis equal;
% hold on;
% 
% flow_graph_between_jump_cones.plot();

%%


% pwintz.plots.namedFigure("Contracted Flow Graphs");
% clf();
% xlim(2.5*[-1, 1]);
% ylim(1.0*[-1, 1]);
% axis equal;
% hold on;
% 
% ctg.plot();
% 

% ctg = union(flow_graph_between_jump_cones, jump_graph)

return

% flow_graph.contractEdgesToBetweenCones(flow_set_1_cone_ndxs, )





% conic_abstraction = ConicAbstraction(conical_partitions, flow_specifications, jump_specifications)


%%
jump_graph_1_mapsto_1 = JumpTransitionGainDigraph(conical_partition_1, jump_set_1_mapsto_1_cone_ndxs, jump_map_matrix_1_mapsto_1);
jump_graph_2_mapsto_2 = JumpTransitionGainDigraph(conical_partition_2, jump_set_2_mapsto_2_cone_ndxs, jump_map_matrix_2_mapsto_2);

pwintz.plots.namedFigure("Jumps Graph 1 to 1");
clf();
xlim(1.2*[-1, 1]);
ylim(1.2*[-1, 1]);
axis square;
hold on;
conical_partition_1.plot();
% conical_partition_1.plotCone(flow_set_cone_ndxs_1, "FaceColor", [0.5, 0.5, 1]);
conical_partition_1.plotCone(jump_graph_1_mapsto_1.jump_set_cone_ndxs,       "FaceColor", [1, 0, 0], "FaceAlpha", 0.5);
conical_partition_1.plotCone(jump_graph_1_mapsto_1.jump_set_image_cone_ndxs, "FaceColor", [1, 0.9, 0.2], "FaceAlpha", 0.5);
jump_graph_1_mapsto_1.plot();
cone_points = jump_graph_1_mapsto_1.getConeRow(jump_set_1_mapsto_1_cone_ndxs).PlotPosition';
pwintz.plots.plotPoints(jump_map_matrix_1_mapsto_1 * cone_points, plotArgs={"Marker", "o", "Color", "red"});

pwintz.plots.namedFigure("Jumps Graph 2 to 2");
clf();
xlim(1.2*[-1, 1]);
ylim(1.2*[-1, 1]);
axis square;
hold on;
conical_partition_2.plot();
% conical_partition_2.plotCone(flow_set_cone_ndxs_2, "FaceColor", [0.5, 0.5, 1]);
conical_partition_2.plotCone(jump_graph_2_mapsto_2.jump_set_cone_ndxs,       "FaceColor", [1, 0, 0], "FaceAlpha", 0.5);
conical_partition_2.plotCone(jump_graph_2_mapsto_2.jump_set_image_cone_ndxs, "FaceColor", [1, 0.9, 0.2], "FaceAlpha", 0.5);
jump_graph_2_mapsto_2.plot();
cone_points = jump_graph_2_mapsto_2.getConeRow(jump_set_2_mapsto_2_cone_ndxs).PlotPosition';
pwintz.plots.plotPoints(jump_map_matrix_2_mapsto_2 * cone_points, plotArgs={"Marker", "o", "Color", "red"});


%%
% flow_transition_gains_digraphs_1 = FlowTransitionGainDigraph(conical_partition_1, flow_set_1_cone_ndxs, flow_map_matrix_1);
% flow_transition_gains_digraphs_2 = FlowTransitionGainDigraph(conical_partition_2, flow_set_2_cone_ndxs, flow_map_matrix_2);


%% 
% conical_partition = ConicalPartition.fromAngles2D(all_angles);

pwintz.plots.namedFigure("Flow Transitions in Mode 1");
clf();
xlim(1.0*[-1, 1]);
ylim(1.0*[-1, 1]);
axis square;
hold on;
conical_partition_1.plot();
conical_partition_1.plotCone(flow_set_1_cone_ndxs, "FaceColor", [0.5, 0.5, 1]);
conical_partition_1.plotCone(jump_set_1_mapsto_1_cone_ndxs, "FaceColor", [1, 0.5, 0.2]);
conical_partition_1.plotCone(jump_set_1_mapsto_2_cone_ndxs, "FaceColor", [1, 0.2, 0.5]);
conical_partition_1.plotCone(jump_graph_1_mapsto_1.jump_set_image_cone_ndxs, "FaceColor", [0.5, 0.5, 1]);
% conical_partition_1.plotCone(jump_graph_2_mapsto_1.jump_set_image_cone_ndxs, "FaceColor", [0.5, 0.5, 1]); % TODO
flow_transition_gains_digraphs_1.plot();
pwintz.plots.plotLinearVectorField(flow_map_matrix_1);
% legend();

%%
pwintz.plots.namedFigure("Flow Transitions in Mode 2");
clf();
xlim(1.0*[-1, 1]);
ylim(1.0*[-1, 1]);
axis square;
hold on;
conical_partition_2.plot();
conical_partition_2.plotCone(flow_set_2_cone_ndxs, "FaceColor", [0.5, 0.5, 1]);
flow_transition_gains_digraphs_2.plot();
pwintz.plots.plotLinearVectorField(flow_map_matrix_1);

%%

warning("Returning early before the end of run_conic_abstraction_with_modes.");
return % TODO: Remove


conic_abstraction_1 = ConicAbstraction(conical_partition_1, flow_set_cone_ndxs, jump_set_cone_ndxs, jump_set_image_cone_ndxs, A_c, A_d);

warning("Returning early before the end of run_conic_abstraction_with_modes.");
return % TODO: Remove

conic_abstraction_1_mapsto_1 = ConicAbstraction.fromAngles2D(...
  flowMapMatrix=-eye(2), ...
  jumpMapMatrix=jump_map_matrix_1_mapsto_1, ...
  jumpSetAngles=jump_set_angles_1_mapsto_1, ...
  flowSetAngles=flow_set_angles_1);
conic_abstraction_1_mapsto_2 = ConicAbstraction.fromAngles2D(...
  flowMapMatrix=-eye(2), ...
  jumpMapMatrix=jump_map_matrix_1_mapsto_2, ...
  jumpSetAngles=jump_set_angles_1_mapsto_2, ...
  flowSetAngles=flow_set_angles_1);
conic_abstraction_2_mapsto_1 = ConicAbstraction.fromAngles2D(...
  flowMapMatrix=-eye(2), ...
  jumpMapMatrix=jump_map_matrix_2_mapsto_1, ...
  jumpSetAngles=jump_set_angles_2_mapsto_1, ...
  flowSetAngles=flow_set_angles_2);
conic_abstraction_2_mapsto_2 = ConicAbstraction.fromAngles2D(...
  flowMapMatrix=-eye(2), ...
  jumpMapMatrix=jump_map_matrix_2_mapsto_2, ...
  jumpSetAngles=jump_set_angles_2_mapsto_2, ...
  flowSetAngles=flow_set_angles_2);

% function 