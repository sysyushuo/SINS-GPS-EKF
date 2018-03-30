%% normalize the quaternion
function quaternion = normalizeQuaternion(quat)

quaternion = zeros(4,1);
length = norm(quat,2);
for i=1:4
    quaternion(i) = quat(i) / length;
end

end