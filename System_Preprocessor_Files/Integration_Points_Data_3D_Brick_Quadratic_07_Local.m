function[W,shape]=Integration_Points_Data_3D_Brick_Quadratic_07_Local(INTP,W)

%% ----------- Region Integration parameters -------------------

        nnp=8;
        ngp=size(INTP,1);
        shape=zeros(ngp,nnp);      
           
    % Generating Shape and dshape matrix
        for i=1:ngp
             x=INTP(i,1);
             y=INTP(i,2);
             z=INTP(i,3);
             shape(i,1:nnp)=Shape_Functions(x,y,z);
        end
                                 
  
        
end
