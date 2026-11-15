function [F_Area, Info] = Ground_Radiation_Precompute(Element_Information, Mesh_Information, FE_Information_Surface, Exterior_Surface_Group, Up_Vector, Dome_Center)


  # Exterior_Surface_Group = 2;
  # Up_Vector = [0 1 0];

  S_Type = Mesh_Information.S_Type;

  F_Area = zeros(Mesh_Information.nnp_max_Model,size(S_Type,2));



  for sN= 1:size(S_Type,2)

      if S_Type(1,sN) == Exterior_Surface_Group

          e = S_Type(2,sN);
          nn = Mesh_Information.Type(2,e);

          gps = FE_Information_Surface.ngp_Surf(sN);

          for gp = 1:gps

              normal = FE_Information_Surface.dTeS_gp(:,gp,sN);
              VF_ground = (1 - dot(normal, Up_Vector))/2;

              T_surf_gp = FE_Information_Surface.T_surf_gp(1:nn,gp,sN);

              F_Area(1:nn,sN) = F_Area(1:nn,sN) + T_surf_gp*VF_ground;

          end

       end

  end
