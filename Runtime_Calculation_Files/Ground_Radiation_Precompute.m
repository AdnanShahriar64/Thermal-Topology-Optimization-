function F_Area=Ground_Radiation_Precompute(Mesh_Information,FE_Information_Surface,Up_Vector)
% Nodal exterior area weighted by view factor to infinite flat ground.
% Exterior group = 1; the dome is centered at the origin.

S_Type=Mesh_Information.S_Type;
up=Up_Vector/norm(Up_Vector);
F_Area=zeros(Mesh_Information.nnp_Model,1);

for sN=1:size(S_Type,2)
    if S_Type(1,sN)==1
        e=S_Type(2,sN);
        nn=Mesh_Information.Type(2,e);
        node_ids=Mesh_Information.IEN(1:nn,e);
        outward=mean(Mesh_Information.Nodes(node_ids,2:4),1)';

        for gp=1:FE_Information_Surface.ngp_Surf(sN)
            normal=FE_Information_Surface.dTeS_gp(:,gp,sN);
            normal=normal/norm(normal);
            if dot(normal,outward)<0
                normal=-normal;
            end
            VF_ground=(1-dot(normal,up))/2;

            F_Area(node_ids)=F_Area(node_ids)+ ...
                FE_Information_Surface.T_surf_gp(1:nn,gp,sN)*VF_ground;
        end
    end
end

end
