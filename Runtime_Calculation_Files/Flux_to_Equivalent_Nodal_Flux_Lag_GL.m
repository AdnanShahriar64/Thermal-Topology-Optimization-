function[QS]=Flux_to_Equivalent_Nodal_Flux_Lag_GL(Element_Information,Mesh_Information,FE_Information_Surface,Nodal_Flux)

%%
    SS=size(Mesh_Information.S_Type);

    ndof=1;

    % SP=size(FE_Information_Surface.P_surf);

    QS=zeros(Mesh_Information.nnp_Model*ndof,1);

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
            
            QV=Nodal_Flux(Mesh_Information.IEN(1:nn,e));

            % QV=zeros(nn,1);
            % PesV(1:3:end-2)=Q_nn;
            
            Qe=FE_Information_Surface.T_surf(1:nn*ndof,count,SN).*QV;                  % Determine Equivalent Nodal FV

            QS(Mesh_Information.IEN(1:nn,e))=QS(Mesh_Information.IEN(1:nn,e))+Qe;     % Put into Total Force Vector

            SNO=SN;

            count=count+1;

    end
   

end
