function f_ground_IR = Ground_Radiation_IR_Load( ...
    F_Area, Wall_Emissivity, Ground_Emissivity, Stefan_Boltzmann, T_ground_IR)
%GROUND_RADIATION_IR_LOAD Ground-origin thermal IR equivalent nodal load [W].
%
%   f_i = epsilon_wall * epsilon_ground * sigma * T_ground_IR^4
%         * integral_Gamma N_i F_ground dGamma
%
%   Assumptions: diffuse-gray ground, uniform ground IR temperature,
%   opaque gray wall with IR absorptivity equal to IR emissivity, no
%   longwave multiple reflections.  T_ground_IR must be in kelvin.

    if ~isvector(F_Area) || any(~isfinite(F_Area))
        error('GroundRadiation:InvalidUnitLoad', 'F_Area must be a finite vector.');
    end
    if ~isscalar(Wall_Emissivity) || ~isfinite(Wall_Emissivity) || ...
            Wall_Emissivity < 0 || Wall_Emissivity > 1 || ...
            ~isscalar(Ground_Emissivity) || ~isfinite(Ground_Emissivity) || ...
            Ground_Emissivity < 0 || Ground_Emissivity > 1
        error('GroundRadiation:InvalidEmissivity', 'Emissivities must be between 0 and 1.');
    end
    if ~isscalar(Stefan_Boltzmann) || ~isfinite(Stefan_Boltzmann) || Stefan_Boltzmann <= 0
        error('GroundRadiation:InvalidSigma', 'Stefan_Boltzmann must be positive.');
    end
    if ~isscalar(T_ground_IR) || ~isfinite(T_ground_IR) || T_ground_IR < 0
        error('GroundRadiation:InvalidGroundTemperature', 'T_ground_IR must be nonnegative kelvin.');
    end

    f_ground_IR = F_Area(:) * ...
        (Wall_Emissivity * Ground_Emissivity * Stefan_Boltzmann * T_ground_IR^4);
end
