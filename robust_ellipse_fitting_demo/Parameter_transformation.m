function [h_hat,k_hat,a_hat,b_hat,theta_hat]=Parameter_transformation (aa_bar)
% % 参数转换
A = aa_bar(1);
B = aa_bar(2);
C = aa_bar(3);
D = aa_bar(4);
E = aa_bar(5);
F = aa_bar(6);
% ============================================================对theta进行值域扩展
% if(aa_bar(1)-aa_bar(3)<0)
%     if(aa_bar(6)>0)
%         theta_hat=(1/2)*atan2(-B,-(A-C));
%     else
%         theta_hat=(1/2)*atan2(B,(A-C));
%     end
% else
%     if(aa_bar(6)>0)
%         theta_hat=(1/2)*atan2(-B,-(A-C));
%     else
%         theta_hat=(1/2)*atan2(B,(A-C));
%     end
% end
    theta_hat=(1/2)*atan2(B,(A-C));
% ============================================================
hk=inv([-2*A -B;-B -2*C])*[D;E];
h_hat=hk(1,:);
k_hat=hk(2,:);
a_hat=sqrt(([h_hat;k_hat]'*[A B/2;B/2  C]*[h_hat;k_hat]-F)/(A*(cos(theta_hat))^2+B*sin(theta_hat)*cos(theta_hat)+C*(sin(theta_hat))^2));
b_hat=sqrt(([h_hat;k_hat]'*[A B/2;B/2  C]*[h_hat;k_hat]-F)/(A*(sin(theta_hat))^2-B*sin(theta_hat)*cos(theta_hat)+C*(cos(theta_hat))^2));
if(abs(a_hat)<abs(b_hat))
    theta_hat=(1/2)*atan2(-B,-(A-C));
    hk=inv([-2*A -B;-B -2*C])*[D;E];
    h_hat=hk(1,:);
    k_hat=hk(2,:);
    a_hat=sqrt(([h_hat;k_hat]'*[A B/2;B/2  C]*[h_hat;k_hat]-F)/(A*(cos(theta_hat))^2+B*sin(theta_hat)*cos(theta_hat)+C*(sin(theta_hat))^2));
    b_hat=sqrt(([h_hat;k_hat]'*[A B/2;B/2  C]*[h_hat;k_hat]-F)/(A*(sin(theta_hat))^2-B*sin(theta_hat)*cos(theta_hat)+C*(cos(theta_hat))^2));
end
end