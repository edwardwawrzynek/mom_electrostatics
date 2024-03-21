% Methods of Moments for Electrostatics
% APPM3310 Final Project
% Edward Wawrzynek, Max Eaton, Andrew Zirger

function mom
    circle = circleMesh(10);
    circle.plotMesh();
    hold on;
    mid = Mesh(circle.midpointSet());
    mid.plotMesh();
end

% construct a mesh for a circle of radius 1 from the specified number of
% points
function mesh = circleMesh(num_pts)
    pts = zeros(num_pts, 2);
    for i = 1:1:num_pts
        pts(i,:) = [cos(i/num_pts * 2*pi) sin(i/num_pts * 2*pi)];
    end

    mesh = Mesh(pts);
end