function [semi_major_axis, eccentricity] = orbital_elements_relative(positions, velocities, masses, G, primary_index)
%ORBITAL_ELEMENTS_RELATIVE Compute osculating a and e relative to a primary.
%
% The calculation treats each secondary body instantaneously as a two-body
% system with the selected primary. This is useful for tracking how the
% N-body perturbations change the moons' osculating orbital elements.

    if nargin < 5
        primary_index = 1;
    end

    n_bodies = size(positions, 1);
    semi_major_axis = NaN(n_bodies, 1);
    eccentricity = NaN(n_bodies, 1);

    for i = 1:n_bodies
        if i == primary_index
            continue;
        end

        r_vec = positions(i, :) - positions(primary_index, :);
        v_vec = velocities(i, :) - velocities(primary_index, :);
        r = norm(r_vec);
        v_sq = dot(v_vec, v_vec);
        mu = G * (masses(primary_index) + masses(i));

        specific_energy = 0.5 * v_sq - mu / r;
        semi_major_axis(i) = -mu / (2 * specific_energy);

        e_vec = ((v_sq - mu / r) * r_vec - dot(r_vec, v_vec) * v_vec) / mu;
        eccentricity(i) = norm(e_vec);
    end
end
