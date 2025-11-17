function dist = createVoxelGrid(mask, stepSize, isocenter)

if ndims(mask)==3
    dist = ones([size(mask), 3]);
else % 2D case
    dist = ones([size(mask), 1, 3]); % singleton for slice
end
if nargin < 3
    isocenter = round(size(dist)/2);
end
if nargin < 2
    stepSize = [1,1,1];
end


for x = 1:size(dist,1)
    dist(x,:,:,1) = (x-isocenter(1)) * stepSize(1);
end
for y = 1:size(dist,2)
    dist(:,y,:,2) = (y-isocenter(2)) * stepSize(2);
end
for z = 1:size(dist,3)
    dist(:,:,z,3) = (z-isocenter(3)) * stepSize(3);
end

end