function [PA, Pb, dPA, dPb] = differentiable_polytope_PH(c, Gdep, Gindep, dc, dGdep, numSliceVar)
% polytope - Converts a zonotope from a G- to a H-representation
% with respect to Z's gradient

indepComb = combinator(size(Gindep,2),2,'c');
    
indepQ = [Gindep(:, indepComb(:, 1)); Gindep(:, indepComb(:, 2))];

% cross product
indepC = [indepQ(2,:).*indepQ(6,:) - indepQ(3,:).*indepQ(5,:); 
          indepQ(3,:).*indepQ(4,:) - indepQ(1,:).*indepQ(6,:); 
          indepQ(1,:).*indepQ(5,:) - indepQ(2,:).*indepQ(4,:)];
normindepC = vecnorm(indepC);
indepC(:,normindepC == 0) = [];
indepC = (indepC ./ normindepC(normindepC > 0))';

depC = [Gdep(2,:).*Gindep(3,:) - Gdep(3,:).*Gindep(2,:); 
        Gdep(3,:).*Gindep(1,:) - Gdep(1,:).*Gindep(3,:); 
        Gdep(1,:).*Gindep(2,:) - Gdep(2,:).*Gindep(1,:)];
    
if nargout > 2
    d_depC = (zeros(size(depC,1), size(depC,2), numSliceVar));
    for i = 1:numSliceVar
        d_depC(:,:,i) = [dGdep(2,i).*Gindep(3,:) - dGdep(3,i).*Gindep(2,:); 
                         dGdep(3,i).*Gindep(1,:) - dGdep(1,i).*Gindep(3,:); 
                         dGdep(1,i).*Gindep(2,:) - dGdep(2,i).*Gindep(1,:)];
    end
end

normdepC = vecnorm(depC);

if nargout > 2
    d_normed_depC = (zeros(size(depC,1), size(depC,2), numSliceVar));

    for i = 1:size(depC,2)
        d_normeddepC_d_depC = [depC(2,i)^2 + depC(3,i)^2, -depC(1,i) * depC(2,i), -depC(1,i) * depC(3,i);
                              -depC(2,i) * depC(1,i), depC(1,i)^2 + depC(3,i)^2, -depC(2,i) * depC(3,i);
                              -depC(3,i) * depC(1,i), -depC(3,i) * depC(2,i), depC(1,i)^2 + depC(2,i)^2] / normdepC(i)^3;

        d_normed_depC(:,i,:) = d_normeddepC_d_depC * squeeze(d_depC(:,i,:));
    end

    d_depC = permute(d_normed_depC, [2,1,3]);
end
                   
depC = (depC./normdepC)';

C = [depC; indepC];

if nargout > 2
    dC = [d_depC; zeros(size(indepC,1), size(indepC,2), numSliceVar)];
end

G = [Gdep, Gindep];

if nargout > 2
    dG = zeros(size(Gindep,1), size(Gindep,2) + 1, numSliceVar);
    dG(:,1,:) = dGdep;
end

deltaD = sum(abs((C*G)'))';

if nargout > 2
    dCG = einsum(C, dG, 'ij,jkl->ikl') + einsum(dC, G, 'ijk,jl->ilk');
    d_deltaD = squeeze(sum(dCG .* sign((C*G)),2));
    d_d = einsum(dC, c, 'ijk,jl->ikl') + C * dc;
end

d = C*c;

PA = [C; -C];
Pb = [d + deltaD; -d + deltaD];

if nargout > 2
    dPA = [dC; -dC];
    dPb = [d_d + d_deltaD; -d_d + d_deltaD];
end
