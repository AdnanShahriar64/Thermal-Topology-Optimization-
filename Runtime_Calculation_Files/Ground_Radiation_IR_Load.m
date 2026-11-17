function f_ground_IR = Ground_Radiation_IR_Load(F_Area, T_ground_IR, Material_Information)

    Wall_Emissivity = Material_Information.epsilon;

    Ground_Emissivity = Material_Information.epsilon_gr;

    Stefan_Boltzmann = Material_Information.StefBC;

    f_ground_IR = F_Area(:) *(Wall_Emissivity * Ground_Emissivity * Stefan_Boltzmann * T_ground_IR^4);
end
