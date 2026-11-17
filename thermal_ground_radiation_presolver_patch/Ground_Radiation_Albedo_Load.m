function f_albedo = Ground_Radiation_Albedo_Load( ...
    F_Area, Solar_Absorptivity, Ground_Albedo, SF, Solar_Ray_Direction, Up_Vector)
%GROUND_RADIATION_ALBEDO_LOAD Reflected-solar equivalent nodal load [W].
%
%   f_i = alpha_wall * rho_ground * SF * max(0,mu_sun)
%         * integral_Gamma N_i F_ground dGamma
%
%   SF is direct-normal solar irradiance [W/m^2].
%   Solar_Ray_Direction is the direction sunlight TRAVELS toward the dome,
%   i.e. the existing solver's Sol_Vec, NOT the vector toward the Sun.
%   mu_sun = -dot(normalized(Solar_Ray_Direction), Up_Vector).
%
%   Assumes a diffuse Lambertian infinite, flat ground, zero shadowing,
%   and uniform albedo; this is an initial approximation.

    if ~isvector(F_Area) || any(~isfinite(F_Area))
        error('GroundRadiation:InvalidUnitLoad', 'F_Area must be a finite vector.');
    end
    if ~isscalar(Solar_Absorptivity) || ~isfinite(Solar_Absorptivity) || ...
            Solar_Absorptivity < 0 || Solar_Absorptivity > 1 || ...
            ~isscalar(Ground_Albedo) || ~isfinite(Ground_Albedo) || ...
            Ground_Albedo < 0 || Ground_Albedo > 1
        error('GroundRadiation:InvalidAlbedo', 'Absorptivity and albedo must be between 0 and 1.');
    end
    if ~isscalar(SF) || ~isfinite(SF) || SF < 0
        error('GroundRadiation:InvalidSolarFlux', 'SF must be nonnegative W/m^2.');
    end
    Solar_Ray_Direction = Solar_Ray_Direction(:);
    Up_Vector = Up_Vector(:);
    if numel(Solar_Ray_Direction) ~= 3 || numel(Up_Vector) ~= 3 || ...
            any(~isfinite(Solar_Ray_Direction)) || any(~isfinite(Up_Vector)) || ...
            norm(Solar_Ray_Direction) == 0 || norm(Up_Vector) == 0
        error('GroundRadiation:InvalidSunDirection', ...
            'Solar_Ray_Direction and Up_Vector must be nonzero finite 3-vectors.');
    end
    sun_ray = Solar_Ray_Direction / norm(Solar_Ray_Direction);
    up = Up_Vector / norm(Up_Vector);

    cos_solar_zenith = max(0, -dot(sun_ray, up));
    f_albedo = F_Area(:) * ...
        (Solar_Absorptivity * Ground_Albedo * SF * cos_solar_zenith);
end
