% classdef ConicalPartitionVertexGraph < handle
%   
%   properties(SetAccess = immutable)
%     % Define private variables.
%     
%     nodes = ones(2, 5);
%     edges = [1, 2]
%     faces = [1, 2, 3]
% 
%     % "vertex_graph" is an undirected vertex_graph of the vertices in the ConicalPartion with edges between if they both border the same cone.
%     vertex_graph; 
%     % cone_graph;   % An undirected vertex_graph of the cones in the ConicalPartion with edges between if they share a vertex (beside zero).
%     vertex_table;
%   end % End of private properties block
% 
%   methods
%     
%     % Constructor
%     function this = ConicalPartitionVertexGraph(vertices)
%       
%     end
%   end
% end % End methods block