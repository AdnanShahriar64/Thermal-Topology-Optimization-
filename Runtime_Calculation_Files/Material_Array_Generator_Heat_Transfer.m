function[Dm]=Material_Array_Generator_Heat_Transfer(Material_Information)

Line='---------------------------------------------------------------------------------------------';
% fprintf('\n%s\n',Line)
% fprintf('%s','6. Calculating material matrix...');

Kx=Material_Information.Kx;

%% Material Matrix initialize
    Dm=zeros(3,3,length(Kx));
    
%% generate material matrix
for i=1:length(Kx)
    

    Dm=[Kx  0   0
        0   Kx  0
        0   0   Kx];
 
 
end

%%

% fprintf('%s\n','Done');
% fprintf('%s\n',Line)

end
