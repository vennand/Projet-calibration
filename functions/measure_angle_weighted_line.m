function data = measure_angle_weighted_line(data)
xyz = data.markers(:,1:4,1)';
x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);

r = mean(xyz);
xyz_fit = bsxfun(@minus,xyz,r);
[~,~,V] = svd(xyz_fit,0);
x_fit = r(1)+(z-r(3))/V(3,1)*V(1,1);
y_fit = r(2)+(z-r(3))/V(3,1)*V(2,1);
% figure(10),clf(10)
% plot3(x,y,z,'b')
% hold on
% plot3(x_fit,y_fit,z,'r')
% scatter3(data.markers(1,:,1), data.markers(2,:,1), data.markers(3,:,1))

xyz_fit = [x_fit, y_fit, z];

weighted_line = xyz_fit(4,:) - xyz_fit(1,:);

angle_deviation = angle_between_vectors(weighted_line, data.gravity');

data.angle_weighted_line = angle_deviation;
end