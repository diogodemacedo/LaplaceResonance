function result = simulate_laplace(europa_mass_factor, n_ganymede_orbits, dt, write_every)
%SIMULATE_LAPLACE Simulate Jupiter, Io, Europa, and Ganymede.
%
% This is a self-contained reconstruction of the numerical model described
% in the 2019 report "Three-Body Laplace Resonance in Jupiter's Moons".
% The original source files are no longer available; the report appendix
% contains the main Octave loop but not the helper functions it called.
%
% Inputs
%   europa_mass_factor : multiplier applied to Europa's physical mass
%   n_ganymede_orbits  : duration in approximate Ganymede orbital periods
%   dt                 : integrator step in normalized time units
%   write_every        : save one sample every N integration steps
%
% Output
%   result : struct containing sampled state and osculating elements

    if nargin < 1, europa_mass_factor = 1; end
    if nargin < 2, n_ganymede_orbits = 40; end
    if nargin < 3, dt = 1.0e-3; end
    if nargin < 4, write_every = 50; end

    % Physical constants and normalization choices used in the report/code.
    G_SI = 6.67408e-11;          % m^3 kg^-1 s^-2
    M0 = 8.9319e22;              % Io mass, kg
    L0 = 4.271e8;                % Io-Jupiter distance scale, m
    T0 = 1.769137786 * 24 * 3600; % Io orbital period, s

    % Convert G into normalized units.
    G = G_SI * M0 * T0^2 / L0^3;

    % Normalized masses.
    m_jupiter  = 1.898e27 / M0;
    m_io       = 8.9319e22 / M0;
    m_europa   = europa_mass_factor * 4.7998e22 / M0;
    m_ganymede = 1.4819e23 / M0;
    masses = [m_jupiter; m_io; m_europa; m_ganymede];

    % Normalized orbital radii used in the report.
    R_io = 1.00;
    R_europa = 1.57;
    R_ganymede = 2.51;

    % Initial positions are coplanar and collinear, reproducing the simple
    % setup described in the report.
    positions = [0, 0, 0;
                 R_io, 0, 0;
                 R_europa, 0, 0;
                 R_ganymede, 0, 0];

    % Circular-orbit speed approximation around Jupiter.
    velocities = zeros(4, 3);
    velocities(2, 2) = sqrt(G * m_jupiter / R_io);
    velocities(3, 2) = sqrt(G * m_jupiter / R_europa);
    velocities(4, 2) = sqrt(G * m_jupiter / R_ganymede);

    % Remove net center-of-mass velocity.
    total_momentum = sum(repmat(masses, 1, size(velocities, 2)) .* velocities, 1);
    v_cm = total_momentum / sum(masses);
    velocities = velocities - repmat(v_cm, size(velocities, 1), 1);

    % Approximate Ganymede period in normalized units and total runtime.
    ganymede_period = 2 * pi * sqrt(R_ganymede^3 / (G * (m_jupiter + m_ganymede)));
    runtime = n_ganymede_orbits * ganymede_period;
    n_steps = ceil(runtime / dt);

    % Small numerical softening; negligible at the orbital scales used here.
    softening = 1.0e-9;

    % Preallocate sampled output.
    max_samples = floor(n_steps / write_every) + 2;
    time = zeros(max_samples, 1);
    sampled_positions = zeros(max_samples, 4, 3);
    semi_major_axis = NaN(max_samples, 4);
    eccentricity = NaN(max_samples, 4);

    acceleration = gravity_acceleration(positions, masses, G, softening);
    [a0, e0] = orbital_elements_relative(positions, velocities, masses, G, 1);

    sample_index = 1;
    time(sample_index) = 0;
    sampled_positions(sample_index, :, :) = reshape(positions, [1, 4, 3]);
    semi_major_axis(sample_index, :) = a0.';
    eccentricity(sample_index, :) = e0.';

    % Velocity-Verlet / leapfrog-equivalent integration.
    for step = 1:n_steps
        velocities = velocities + 0.5 * acceleration * dt;
        positions = positions + velocities * dt;

        new_acceleration = gravity_acceleration(positions, masses, G, softening);
        velocities = velocities + 0.5 * new_acceleration * dt;
        acceleration = new_acceleration;

        if mod(step, write_every) == 0 || step == n_steps
            sample_index = sample_index + 1;
            [a_now, e_now] = orbital_elements_relative(positions, velocities, masses, G, 1);

            time(sample_index) = min(step * dt, runtime);
            sampled_positions(sample_index, :, :) = reshape(positions, [1, 4, 3]);
            semi_major_axis(sample_index, :) = a_now.';
            eccentricity(sample_index, :) = e_now.';
        end
    end

    % Trim unused preallocated rows.
    time = time(1:sample_index);
    sampled_positions = sampled_positions(1:sample_index, :, :);
    semi_major_axis = semi_major_axis(1:sample_index, :);
    eccentricity = eccentricity(1:sample_index, :);

    result = struct();
    result.europa_mass_factor = europa_mass_factor;
    result.G = G;
    result.mass_scale_kg = M0;
    result.length_scale_m = L0;
    result.time_scale_s = T0;
    result.ganymede_period = ganymede_period;
    result.runtime = runtime;
    result.dt = dt;
    result.masses = masses;
    result.time = time;
    result.positions = sampled_positions;
    result.semi_major_axis = semi_major_axis;
    result.eccentricity = eccentricity;
end
