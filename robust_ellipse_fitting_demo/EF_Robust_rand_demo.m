clc;clear;close all;
z=50; %Number of random scenes
It=50; %Number of iterations
V=500; %Number of Monte Carlo runs
N_n=100; %Number of data points
num=0;
tic
for outlier=40 %Number of outliers
    num=num+1
%%==================================Generating ellipse parameters
rand('seed',10);
h_all = 20*rand(z,1);
k_all = 20*rand(z,1);
b_all = 20*(2*rand(z,1)-ones(z,1))+30;
sigma_all = 0.005*b_all;
theta_all=90*(2*rand(z,1)-ones(z,1))/360*2*pi;
alpha_true_all = cell(z,1);
M_all= cell(z,1);
M1_all= cell(z,1);
MM_all= cell(z,1);
MM1_all= cell(z,1);
randn('seed',160);
for i = 1:z
a_all(i,1) =((50-b_all(i))/2)*(2*rand(1,1)-ones(1,1))+((60+b_all(i))/2);
alpha_true_all{i} = rand(N_n,1)*2*pi;
M_all{i} = randn(N_n,V);
M1_all{i} = randn(N_n,V);
MM_all{i} = b_all(i)*(2*rand(outlier,V)-ones(outlier,V));
MM1_all{i} = b_all(i)*(2*rand(outlier,V)-ones(outlier,V));
end
%%==================================
for ii=50 %%Scene 50
sigma=sigma_all(ii);
h = h_all(ii);
k = k_all(ii);
b = b_all(ii);
a = a_all(ii);
theta=theta_all(ii);
noise = N_n-outlier;
alpha_true = alpha_true_all{ii};
M = M_all{ii};
M1 = M1_all{ii};


MM = MM_all{ii};
MM1 = MM1_all{ii};
alpha_true0 = 0:2*pi/2000:2*pi;
Np = length(alpha_true0);
for i=1:Np
    cos_alpha(i) = sqrt(1)*cos(alpha_true0(i));
    sin_alpha(i) = sqrt(1)*sin(alpha_true0(i));
x(i) = (a*cos_alpha(i)*cos(theta) - b*sin_alpha(i)*sin(theta)) + h;
y(i) = (a*cos_alpha(i)*sin(theta) + b*sin_alpha(i)*cos(theta)) + k;
end
% plot(x,y,'k'); hold on;
%  =========================================== Generating data points
for v=379 %%The 379th Monte Carlo run
    cos_alpha=[];
    sin_alpha=[];
    n=[];
    xx=[];yy=[];
    for i=1:noise
        ns = sigma*M(i,v);%%Noise error on the xlabel
        ns1 = sigma*M1(i,v);%%Noise error on the ylabel
        cos_alpha(i) = sqrt(1)*cos(alpha_true(i,1));
        sin_alpha(i) = sqrt(1)*sin(alpha_true(i,1));
    xx(i) = (a*cos_alpha(i)*cos(theta) - b*sin_alpha(i)*sin(theta)) + h+ns;
    yy(i) = (a*cos_alpha(i)*sin(theta) + b*sin_alpha(i)*cos(theta)) + k+ns1;
    end
  % % ========================================== outliers
 for i=noise+1:N_n
         ns = sigma*M(i,v);
        ns1 = sigma*M1(i,v);
         NS = MM(i-noise,v); %%Outliers error on the xlabel
        NS1 = MM1(i-noise,v); %%Outliers error on the ylabel
        cos_alpha(i) = sqrt(1)*cos(alpha_true(i,1));
        sin_alpha(i) = sqrt(1)*sin(alpha_true(i,1));
    xx(i) = (a*cos_alpha(i)*cos(theta) - b*sin_alpha(i)*sin(theta)) + h+ns+NS;
    yy(i) = (a*cos_alpha(i)*sin(theta) + b*sin_alpha(i)*cos(theta)) + k+ns1+NS1;
 end     
%     plot(xx,yy,'r+');
for i=noise+1:N_n
    text(xx(i)+0.01,yy(i)+0.01,num2str(i));
end
%  ===========================================
    for i = 1:N_n
        n(:,i)= [xx(i)^2 xx(i)*yy(i) yy(i)^2 xx(i) yy(i) 1].';
    end  
    tic
        [aa,itLPit,ki]=EF_LP_Roubst(n,N_n,xx,yy,It);
        [h_hat1,k_hat1,a_hat1,b_hat1,theta_hat1]=Parameter_transformation (aa);
        LPRMSE=(h_hat1-h)^2+(k_hat1-k)^2+(a_hat1-a)^2+(b_hat1-b)^2+(theta_hat1-theta)^2
    timeLP(num,v,ii)=toc;
end
end
end
toc
% save('rand_cluster_position.mat')
% %=========================================================  fit ellipse 
alpha_hat = 0:2*pi/2000:2*pi;
Np = length(alpha_hat);
for i=1:Np
x_hat1(i) = (a_hat1*cos(alpha_hat(i))*cos(theta_hat1) - b_hat1*sin(alpha_hat(i))*sin(theta_hat1)) + h_hat1;
y_hat1(i) = (a_hat1*cos(alpha_hat(i))*sin(theta_hat1) + b_hat1*sin(alpha_hat(i))*cos(theta_hat1)) + k_hat1;
end

for i=1:noise
    x_edge(i)=xx(i);
    y_edge(i)=yy(i);
end
for i=noise+1:N_n
    x_outlier(i-noise)=xx(i);
    y_outlier(i-noise)=yy(i);
end
figure(3)
plot(x_edge,y_edge,'r+','LineWidth',1.5);hold on;
plot(x_outlier,y_outlier,'ro','LineWidth',1.5);hold on;
plot(x,y,'k-','LineWidth',1.5); hold on;
plot(x_hat1,y_hat1,'m-','LineWidth',1.5);hold on;
% legend('Edge Points','Outliers','True Ellipse','MCC-GSSOCP','MCC-Laplacian','MCC-Gaussian','ADMM','RANSAC','SAREfit','Munoz');
legend('Edge Points','Outliers','True Ellipse','MCC-Laplacian')
xlabel('xlabel');
ylabel('ylabel');
axis([-80 120 -100 100])
