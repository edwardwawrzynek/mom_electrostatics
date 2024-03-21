classdef Mesh
    properties
        num_pts
        points  % array of ordered points in the mesh
        charge_density
    end
    methods
        function obj = Mesh(pts)
            obj.points = pts;
            obj.num_pts = length(pts);
        end
        
        % Display the mesh points and edges
        function plotMesh(obj)
            plot([obj.points(:,1); obj.points(1,1)], [obj.points(:,2); obj.points(1,2)], '-o');
        end
        
        % Generate a set of points that are at the midpoint of the mesh
        % points
        function midpoints = midpointSet(obj)
            midpoints = zeros(obj.num_pts, 2);
            for i = 1:1:obj.num_pts
                midpoints(i,:) = 0.5 * (obj.points(i,:) + obj.points(mod(i, obj.num_pts)+1,:));
            end
        end

    end
end