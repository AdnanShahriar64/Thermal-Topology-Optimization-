% function [F_Area] = Ground_Radiation_Precompute(Mesh_Information, FE_Information_Surface, Exterior_Surface_Group, Up_Vector)


  Exterior_Surface_Group = 1;
  Up_Vector = [0 1 0];

  S_Type = Mesh_Information.S_Type;

  % F_Area = zeros(Mesh_Information.nnp_max_Model,size(S_Type,2));

  % Initialize before the surface loop
  F_Area = zeros(Mesh_Information.nnp_Model,1);



  for sN= 1:size(S_Type,2)

      if S_Type(1,sN) == Exterior_Surface_Group

          e = S_Type(2,sN);
          nn = Mesh_Information.Type(2,e);

          gps = FE_Information_Surface.ngp_Surf(sN);

          node_ids = Mesh_Information.IEN(1:nn,e);

          for gp = 1:gps

              up = Up_Vector(:) / norm(Up_Vector);

              dTe = FE_Information_Surface.dTeS_gp(:,gp,sN);
              normal = dTe / norm(dTe);

              VF_ground = (1 - dot(normal, up)) / 2;

              T_surf_gp = FE_Information_Surface.T_surf_gp(1:nn,gp,sN);

              % F_Area(1:nn,sN) = F_Area(1:nn,sN) + T_surf_gp*VF_ground;

              % Inside the selected face and Gauss-point loops
              F_Area(node_ids) = F_Area(node_ids) + T_surf_gp * VF_ground;

          end

       end

  end
