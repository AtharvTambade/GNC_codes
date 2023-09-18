function quaternion_noise_norm = noise(quaternion,sigma_theta,sigma_vector)
if quaternion(1)==1
    quaternion_noise_norm=[1,0,0,0];
else
    theta=2*acosd(quaternion(1));
    theta_radian=theta*pi/180;
    noise_theta=theta + normrnd(0,sigma_theta);
    factor=sin(theta_radian/2);
    beta_x=acosd(quaternion(2)/factor);
    noise_beta_x=beta_x + normrnd(0,sigma_vector);
    noise_beta_x_radian=noise_beta_x*pi/180;
    beta_y=acosd(quaternion(3)/factor);
    noise_beta_y=beta_y + normrnd(0,sigma_vector);
    noise_beta_y_radian=noise_beta_y*pi/180;
    beta_z=acosd(quaternion(4)/factor);
    noise_beta_z=beta_z + normrnd(0,sigma_vector);
    noise_beta_z_radian=noise_beta_z*pi/180;
    noise_theta_radian=noise_theta*pi/180;
    quaternion_noise=[cos(noise_theta_radian/2),sin(noise_theta_radian/2)*cos(noise_beta_x_radian),sin(noise_theta_radian/2)*cos(noise_beta_y_radian),sin(noise_theta_radian/2)*cos(noise_beta_z_radian)];
    quaternion_noise_norm = quaternion_noise/norm(quaternion_noise);
end
end
    
