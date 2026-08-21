function acceleration = gravity_acceleration(positions, masses, G, softening)
%GRAVITY_ACCELERATION Compute Newtonian N-body accelerations.
%
%   acceleration = gravity_acceleration(positions, masses, G, softening)
%
% Inputs
%   positions : N x D matrix of body positions
%   masses    : N x 1 vector of body masses
%   G         : gravitational constant in the chosen normalized units
%   softening : small length used only to avoid numerical singularities
%
% Output
%   acceleration : N x D matrix of accelerations

    n_bodies = size(positions, 1);
    acceleration = zeros(size(positions));

    for i = 1:n_bodies-1
        for j = i+1:n_bodies
            r_ij = positions(j, :) - positions(i, :);
            distance_sq = sum(r_ij .^ 2) + softening^2;
            inv_distance_cubed = distance_sq^(-1.5);

            acceleration(i, :) = acceleration(i, :) ...
                + G * masses(j) * r_ij * inv_distance_cubed;
            acceleration(j, :) = acceleration(j, :) ...
                - G * masses(i) * r_ij * inv_distance_cubed;
        end
    end
end
