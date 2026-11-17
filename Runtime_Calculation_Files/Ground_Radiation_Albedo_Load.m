function f_albedo = Ground_Radiation_Albedo_Load(F_Area, SF, Solar_Ray_Direction, Up_Vector,Material_Information)

    Solar_Absorptivity = Material_Information.absorbability;

    Ground_Albedo = Material_Information.ground_albedo;


    sun_ray = Solar_Ray_Direction / norm(Solar_Ray_Direction);
    up = Up_Vector / norm(Up_Vector);

    cos_solar_zenith = max(0, -dot(sun_ray, up));
    f_albedo = F_Area(:) * (Solar_Absorptivity * Ground_Albedo * SF * cos_solar_zenith);
end
