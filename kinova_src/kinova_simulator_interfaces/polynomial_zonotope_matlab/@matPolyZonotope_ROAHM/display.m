function display(matpZ)

%------------- BEGIN CODE --------------

%display dimension, generators
disp('dimension: ');
disp(size(matpZ.C));
disp('nr of dep. generators: ');
disp(size(matpZ.G, 3));
disp('nr of ind. generators: ');
disp(size(matpZ.Grest, 3));
%display center
disp('center: ');
disp(matpZ.C);

%display generators
disp('dep. generators: ');
for i=1:size(matpZ.G, 3)
    disp(matpZ.G(:, :, i)); 
    disp('---------------'); 
end

disp('ind. generators: ');
for i=1:size(matpZ.Grest, 3)
    disp(matpZ.Grest(:, :, i)); 
    disp('---------------'); 
end

% display exponential matrix
disp('exponential matrix: ');
disp(matpZ.expMat);

% display id
disp('id:');
disp(matpZ.id);

%------------- END OF CODE --------------