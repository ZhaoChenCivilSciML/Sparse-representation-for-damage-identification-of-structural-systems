function [shapeFEM2,eigvalFEM2,shapeMea2] = MACpair(shapeFEM,eigvalFEM,shapeMea)
% This function
% 1. sorts orders of FEM&Mea mode shapes
% 2. adjusts the orientations of measured mode shapes
% so that mode shapes from FEM&Mea are similar to the largest extent.

% Select and sort the order
for i = 1:size(shapeFEM,2)
    for j = 1:size(shapeMea,2)
        MAC(i,j) = (shapeFEM(:,i)'*shapeMea(:,j))^2/(norm(shapeFEM(:,i)))^2/(norm(shapeMea(:,j)))^2;
    end
end
for i = 1:size(shapeMea,2)
    [~,ind] = max(MAC(:,i));
    shapeFEM2(:,i) = shapeFEM(:,ind);
    eigvalFEM2(i,1) = eigvalFEM(ind);
end
if find(eigvalFEM2 == 0) ~= 0
    error('error ocurrs in MAC!')
end

% Adjust the orientation
shapeMea2 = shapeMea;
for i = 1:size(shapeFEM2,2)
    if norm(shapeFEM2(:,i)-shapeMea(:,i))>norm(-shapeFEM2(:,i)-shapeMea(:,i))
        shapeMea2(:,i) = -shapeMea(:,i);
    end
end