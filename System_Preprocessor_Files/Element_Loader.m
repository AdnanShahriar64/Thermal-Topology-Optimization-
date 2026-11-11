function[Element_Information,Nodes_Elemental,nnp_Elemental,ngp_Elemental,Intp_Elemental,Shape,DShape,W,Surf_Nodes_Elemental,Surf_Intp_Elemental,Surf_ngp_Elemental,Surf_Shape,Surf_Dshape,Surf_W,Par_Coord]=Element_Loader
% 1 - Linear Tetra                  - 4 Node
% 2 - Quadratic Tetra               - 10 Node
% 3 - Linear Wedge                  - 6 Node    - 6  GP
% 4 - Quadratic Wedge               - 15 Node   - 21 GP
% 5 - Pyramid                       - 18 Node   - 21 GP
% 6 - Linear Brick                  - 8 Node    - 8  GP
% 7 - Quadratic Brick               - 20 Node   - 27 GP
% 8 - Enriched Quadratic Brick      - 27 Node   - 27 GP  
% 9 - Manual Element 9 Node Brick   - 9 Node    - 8  GP

%% Element Library

MaxGp=27;
MaxN=20;
MaxS=6;
MaxSGp=4;
Maxnnp=27;

Lib_Dir=fullfile(fileparts(fileparts(mfilename('fullpath'))),'Element_Library');  % <main folder>/Element_Library
list=dir(fullfile(Lib_Dir,'*.m'));

ElL=length(list);

Nodes_Elemental=zeros(Maxnnp,3,ElL);
Intp_Elemental=zeros(MaxSGp,3,ElL);

W=zeros(MaxGp,ElL);
Shape=zeros(MaxGp,MaxN,ElL);
DShape=zeros(3,MaxN,ElL);
nnp_Elemental=zeros(1,ElL);
ngp_Elemental=zeros(1,ElL);

Surf_ngp_Elemental=zeros(6,ElL);
Surf_nnp_Elemental=zeros(6,ElL);
Par_Coord=zeros(3,2,8,ElL);
Surf_W=zeros(MaxSGp,ElL);
Surf_Shape=zeros(MaxSGp,MaxN,MaxS,ElL);
Surf_Dshape=zeros(3,MaxN,MaxSGp,MaxS,ElL);
% Surf_Dshape=zeros(3,1,1,1,ElL);

Surf_Nodes_Elemental=zeros(round(Maxnnp/3),MaxS,ElL);
Surf_Intp_Elemental=zeros(MaxSGp,3,6,ElL);



Line='------------------------------------------------------------------------------------------------------------------------------------------------------------';
fprintf('\n%s\n',Line)
fprintf('Loading elements:\n')
fprintf('%s\n',Line)
fprintf('%60s   %3s     %16s     %19s       %12s       %12s\n','Imported Function Name','DOF','Nodes in Element','Surfaces in Element','Integration points','Element Type') 
fprintf('%s\n',Line)


for i=1:ElL
   
    fcnName=list(i).name;
    fcnName=erase(fcnName,'.m');
    
    %% For Standard FEM
    if fcnName(end)~='X'  
        
        [NodesO,INTPO,WO,ShapeO,DshapeO,Surf_Nodes_ElementalO,Surf_INTPO,Surf_WO,Surf_ShapeO,Surf_DshapeO,Par_CoordO]=eval(fcnName);

        ElNum=str2num(fcnName(end-1:end));
        if size(ElNum)==0
            ElNum=str2num(fcnName(end));
        end

        % Region properties
        SE=size(ShapeO);
        SD=size(DshapeO);
        ngp=SE(1);
        nnp=SE(2);
        ndof=SD(1);
        % ngpS=size(Surf_ShapeO,1);

        W(1:ngp,ElNum)=WO;
        Shape(1:ngp,1:nnp,ElNum)=ShapeO;
        DShape(1:ndof,1:nnp,1:ngp,ElNum)=DshapeO;
        nnp_Elemental(1,ElNum)=nnp;
        ngp_Elemental(1,ElNum)=ngp;
        Nodes_Elemental(1:nnp,:,ElNum)=NodesO;

        Intp_Elemental(1:ngp,:,ElNum)=INTPO;
        
        % Surface properties
%         L=size(Surf_DshapeO);
%         if length(L)<3
%             sn=1;
%             Surf_ngp_ElementalO=1;
%         elseif length(L)<4
%             sn=1;
%         else
%             sn=L(4);
%             Surf_ngp_ElementalO=L(3);
%         end

        sn=size(Surf_DshapeO,4); % Number of surfaces
        S=size(Surf_WO);
        for ss=1:S(2)
            gpsCount=0;
            for gps=1:S(1)
                if Surf_WO(gps,ss)>0
                   gpsCount=gpsCount+1;
                   Surf_W(gpsCount,ss,ElNum)=Surf_WO(gpsCount,ss);
                   Surf_Intp_Elemental(gpsCount,:,ss,ElNum)=Surf_INTPO(gpsCount,:,ss);
                end
                Surf_ngp_Elemental(ss,ElNum)=gpsCount;
            end
        end
        maxSgp=max(Surf_ngp_Elemental(:,ElNum));

        %Surf_Nodes_ElementalO,Surf_INTPO

        nnpS=size(Surf_Nodes_ElementalO,1);

        for ss=1:S(2)
            count=0;
            for nn=1:nnpS
                if Surf_Nodes_ElementalO(nn,ss)>0
                    Surf_nnp_Elemental(ss,ElNum)=count+1;
                    count=count+1;
                end
                
            end
        end
        Surf_Nodes_Elemental(1:nnpS,1:sn,ElNum)=Surf_Nodes_ElementalO;




        Surf_Dshape(1:ndof,1:nnp,1:maxSgp,1:sn,ElNum)=Surf_DshapeO;
        Surf_Shape(1:maxSgp,1:nnp,1:sn,ElNum)=Surf_ShapeO;
        
        Par_Coord(:,:,1:sn,ElNum)=Par_CoordO;

        fprintf('%60s   %3d     %11d     %17d       %15d         %15d\n',fcnName,ndof,nnp,sn,ngp,ElNum)

        % Clearing var
        clear WO shapeO DshapeO Surf_WO Surf_ngp_ElementalO Surf_ShapeO Par_CoordO Surf_DshapeO
    
    end
    
    
end


% Element_Information=struct('Nodes_Elemental','nnp_Elemental','ngp_Elemental','Intp_Elemental','Shape,DShape','W,Surf_Nodes_Elemental','Surf_Intp_Elemental','Surf_ngp_Elemental','Surf_Shape','Surf_Dshape','Surf_W','Par_Coord',[]);

Element_Information.Nodes_Elemental=Nodes_Elemental;
Element_Information.nnp_Elemental=nnp_Elemental;
Element_Information.ngp_Elemental=ngp_Elemental;
Element_Information.Intp_Elemental=Intp_Elemental;
Element_Information.Shape_Elemental=Shape;
Element_Information.DShape_Elemental=DShape;
Element_Information.W_Elemental=W;
Element_Information.Surf_Nodes_Elemental=Surf_Nodes_Elemental;
Element_Information.Surf_Intp_Elemental=Surf_Intp_Elemental;
Element_Information.Surf_ngp_Elemental=Surf_ngp_Elemental;
Element_Information.Surf_nnp_Elemental=Surf_nnp_Elemental;
Element_Information.Surf_Shape_Elemental=Surf_Shape;
Element_Information.Surf_Dshape_Elemental=Surf_Dshape;
Element_Information.Surf_W_Elemental=Surf_W;
Element_Information.Par_Coord_Elemental=Par_Coord;

fprintf('%s\n\n',Line)

end
