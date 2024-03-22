classdef Mesh
    properties
        num_pts
        points  % array of ordered points in the mesh
        basis % the type of basis and testing functions to use
        weights % solved weight of charge distribution basis
    end
    methods
        function obj = Mesh(pts, basis)
            obj.points = pts;
            obj.num_pts = length(pts);
            obj.basis = basis;
            obj.weights = zeros(obj.num_pts,1);
        end
        
        % Display the mesh points and edges
        function plotMesh(obj)
            plot([obj.points(:,1); obj.points(1,1)], [obj.points(:,2); obj.points(1,2)], '-o');
        end

        % Display the mesh colored by charge density
        function plotCharge(obj)
            % compute the absolue maximum charge value
            max_charge = max(obj.weights);

            cm = colormap(jet);
            colorbar;

            plot_pts = [[obj.points(:,1); obj.points(1,1)], [obj.points(:,2); obj.points(1,2)]];
            for i = 1:obj.num_pts
                color0 = cm(floor(0.5*(obj.weights(i)/max_charge + 1.0) * (size(cm, 1)-1)) + 1,:);
                
                line([plot_pts(i,1) plot_pts(i+1,1)], [plot_pts(i,2) plot_pts(i+1,2)], ...
                    'Color', color0, ...
                    'LineWidth', 2);
            end
        end

        % Compute voltage over space created by the charge distribution
        function [v, xleft, xright, yleft, yright] = computeVoltage(obj, scale_factor, n)
            % Get range over which to evaluate voltage
            x = obj.points(:,1);
            y = obj.points(:,2);

            xcenter = 0.5 * (max(x)+min(x));
            xdist = max(x)-min(x);
            ycenter = 0.5 * (max(y)+min(y));
            ydist = max(y)-min(y);

            xleft = xcenter - 0.5 * scale_factor * xdist;
            xright = xcenter + 0.5 * scale_factor * xdist;
            yleft = ycenter - 0.5 * scale_factor * ydist;
            yright = ycenter + 0.5 * scale_factor * ydist;

            % evaluate voltage
            pts_neighbors = [obj.points(obj.num_pts,:); obj.points; obj.points(1,:)];

            v = zeros(n);
            
            % minimum distance at which to evaluate v
            min_dist = 0.5 * min(xdist, ydist)*scale_factor/n;
            
            for i = 1:n
                for j = 1:n
                    pt = [xleft + (i-1) * scale_factor * xdist/n, yleft + (j-1) * scale_factor * ydist/n];
                    % Sum contribution from each basis
                    for k = 1:obj.num_pts
                        %v(i,j) = v(i,j) + obj.weights(k) * obj.basis.evaluateVoltage(pt, pts_neighbors(k+1,:), pts_neighbors(k,:), pts_neighbors(k+2,:), );
                        dist = max(norm(pt - obj.points(k,:)), min_dist);

                        v(i,j) = v(i,j) + 1/(4*pi) * 1 / dist;
                    end
                end
            end
        end

        % Display voltage and electric field
        function plotVoltage(obj, scale_factor, n)
            [v, xleft, xright, yleft, yright] = obj.computeVoltage(scale_factor, n);
            
            figure;
            colormap("jet");
            imagesc(rot90(v));
            title("Voltage");
            colorbar;

            % Compute electric field strength
            v_upper = v(1:end-1,1:end-1);
            E_x = (v(2:end, 1:end-1) - v_upper) / ((xright-xleft)/n);
            E_y = (v(1:end-1, 2:end) - v_upper) / ((yright-yleft)/n);

            E = sqrt(E_x .^ 2 + E_y .^ 2);
            
            figure;
            colormap("jet");
            imagesc(rot90(E));
            colorbar;
            title("Electric Field Intensity [dB]");
        end

        % Solve for the charge distribution created by charging the mesh to
        % V0.
        function obj = solve(obj, V0)
            % Coefficient matrix
            A = zeros(obj.num_pts);
            % weighted voltage vector
            v = zeros(obj.num_pts,1);

            % compute points with neighbors included at beginning and end
            pts_neighbors = [obj.points(obj.num_pts,:); obj.points; obj.points(1,:)];
            
            % compute voltage vector
            for i = 1:1:obj.num_pts
                v(i) = V0 * obj.basis.magnitude(pts_neighbors(i+1,:), pts_neighbors(i, :), pts_neighbors(i+2, :));
            end

            % compute coefficient matrix
            for i = 1:1:obj.num_pts
                for j = 1:1:obj.num_pts
                    A(i,j) = obj.basis.innerProduct(pts_neighbors(i+1, :), pts_neighbors(i, :), pts_neighbors(i+2, :), pts_neighbors(j+1, :), pts_neighbors(j, :), pts_neighbors(j+2, :));
                end
            end
            
            % TODO: investigate how to solve when A is singular
            obj.weights = lsqr(A,v,1e-6,128);
        end
        
        % Generate a set of points that are at the midpoint of the mesh
        % points
        % function midpoints = midpointSet(obj)
        %    midpoints = zeros(obj.num_pts, 2);
        %    for i = 1:1:obj.num_pts
        %        midpoints(i,:) = 0.5 * (obj.points(i,:) + obj.points(mod(i, obj.num_pts)+1,:));
        %    end
        % end

    end
end