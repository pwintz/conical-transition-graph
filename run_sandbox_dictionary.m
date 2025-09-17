
min_gains = [];
max_gains = [];
for gamma = gammas
  % gamma
  % gamma_results{gamma} 

  min_gains(end+1) = gamma_results{gamma}.min_gain;
  max_gains(end+1) = gamma_results{gamma}.max_gain;
end

figure(1);
clf();

plot(gammas, min_gains, "Marker", ".", "MarkerSize", 7, "DisplayName",  "Maximum Cycle Weight");
hold on;
plot(gammas, max_gains, "Marker", ".", "MarkerSize", 7, "DisplayName",  "Minimum Cycle Weight");
plot(gammas([1, end]), [1, 1], '--k', "DisplayName", "Threshold");
lgd = legend();
lgd.Location = "best";
xlabel("\gamma");
ylabel("Weight");
title("Cycle Walk Weights vs. \gamma Parameter");

set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');

pwintz.plots.saveExampleFigure('C:\Users\pwintz\code\Zeno\latex\images', 'numerical_example_weights_vs_gamma', width=500, height=300);



% 
% pwintz.math.angle2UnitVector(pi     )
% pwintz.math.angle2UnitVector(pi + pi/12)
% pwintz.math.angle2UnitVector(pi/2      )
% pwintz.math.angle2UnitVector(2*pi/3)
% pwintz.math.angle2UnitVector(pi - pi/12 )
% pwintz.math.angle2UnitVector(pi)
% pwintz.math.angle2UnitVector(0)
% pwintz.math.angle2UnitVector(pi/12)



% pwintz.math.atan2([-4; -1]) - ( pi + pi/12) + 2*pi
% pwintz.math.atan2([0; 1])- pi/2
% pwintz.math.atan2([-1; 2])- 2*pi/3
% pwintz.math.atan2([-4; 1]) - (pi - pi/12)
% pwintz.math.atan2([-1; 0]) - pi
% pwintz.math.atan2([1; 0]) - 0
% pwintz.math.atan2([4; 1]) - pi/12