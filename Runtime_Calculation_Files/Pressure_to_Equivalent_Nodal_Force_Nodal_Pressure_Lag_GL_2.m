function[fS]=Pressure_to_Equivalent_Nodal_Force_Nodal_Pressure_Lag_GL_2(Element_Information,Mesh_Information,FE_Information_Surface,Nodal_Pressure)

%%
    SS=size(Mesh_Information.S_Type);

    ndof=Mesh_Information.ndof_Model;

    % SP=size(FE_Information_Surface.P_surf);

    fS=zeros(Mesh_Information.nnp_Model*ndof,1);

    count=1;

    SNO=1;
    for i=1:SS(2) % For each surfaces


            SN=Mesh_Information.S_Type(1,i);                            % Extract Surface no 
            if SN>SNO
                count=1;
            end

            Surf_Type=Mesh_Information.S_Type(3,i);
        
            e=Mesh_Information.S_Type(2,i);                         % Extract element no

            nn=Element_Information.nnp_Elemental(Mesh_Information.Type(1,e));

            ng=Element_Information.ngp_Elemental(Mesh_Information.Type(1,e));

            if nn~=ng || max(max(Element_Information.Surf_Shape_Elemental(:,:,Surf_Type,Mesh_Information.Type(1,e))))~=1
                error('Change the Pressure calculation file, this doens seem like a Lagrangian element with GL nodes')
            end

            % Pes=Ps(i);                                       % Extract Element Pressure

            Pes=Nodal_Pressure(Mesh_Information.IEN(1:nn,e));
            
            PesV=zeros(nn*3,1);
            PesV(1:3:end-2)=Pes;
            PesV(2:3:end-1)=Pes;
            PesV(3:3:end-0)=Pes;
            
            fe=FE_Information_Surface.P_surf(1:nn*ndof,count,SN).*PesV;                  % Determine Equivalent Nodal FV

            fS(Mesh_Information.LM_Model(1:nn*3,e))=fS(Mesh_Information.LM_Model(1:nn*3,e))+fe;     % Put into Total Force Vector

            SNO=SN;

            count=count+1;

    end
   

end
