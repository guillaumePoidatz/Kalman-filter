% Model for a fast implementation of Kalman filter
% Ph. Bonnifait
clear
clc
load simulated_data.mat; %load variables in the workspace

nepoch=length(t);

%size of the state
n=3; %for instance X=(x,y,vx,vy) - to be changed

% Variables used to store the results
Xs=zeros(n,nepoch); 
%If we want to store some terms of the covariance matrix, wr allocate
%before some memory
Px1=zeros(1,nepoch); %variance
Px1x2=zeros(1,nepoch); %covariance

% Initial state of the filter (to be changed)
X=zeros(n,1);
P=100*eye(n,n);
Qv = 100000000000000; % standard deviation of the noise on the model (speed)
Qd = 10000; % standard deviation of the noise on the model (distance)
for i=1:nepoch
    % first we compute the matrices of the model
    Q = [dt(i)^2/4*Qv,0,0;0,Qv,0;0,0,Qd]; % random noise on the model
    A = [1,dt(i),0;0,1,0;0,0,1]; % no command
    C = [(X(1)-gps.s1.x(i))/sqrt((gps.s1.x(i)-X(1))^2+gps.s1.y(i)^2),0,1;...
        (X(1)-gps.s2.x(i))/sqrt((gps.s2.x(i)-X(1))^2+gps.s2.y(i)^2),0,1;...
        (X(1)-gps.s3.x(i))/sqrt((gps.s3.x(i)-X(1))^2+gps.s3.y(i)^2),0,1;...
        0,1,0]
    G = [sqrt((gps.s1.x(i)-X(1))^2+gps.s1.y(i)^2)+X(3);sqrt((X(1)-gps.s2.x(i))^2+gps.s2.y(i)^2)+X(3);...
        sqrt((X(1)-gps.s3.x(i))^2+gps.s3.y(i)^2)+X(3);X(2)]
    R = [gps.s1.sPR,0,0,0;0,gps.s2.sPR,0,0;0,0,gps.s3.sPR,0;0,0,0,tachy.sv].^2; % random noise on the measurements
    Y = [gps.s1.PR(i);gps.s2.PR(i);gps.s3.PR(i);tachy.v(i)]
    % measurement update (estimation)
    % error on the estimation
    epsilon = Y-G;
    K = P*C'*1/(C*P*C'+R); % Kalman gain (Pxy)*Py-1
    X = X+K*epsilon; % update of X (Xk|k)
    P = (eye(3)-K*C)*P*(1-K*C)'+K*R*K'; % cov matrix of the update (Pk|k)
    
    % Storage (corresponds to the output of the filter)
    Xs(:,i)=X;
    Px1(i)=P(1,1);
    Px1x2(i)=P(1,2);
    %...
    
    % time update (prediction)
    X = A*X; % prediction (Xk|k-1)

    % covariance matrix of the estimate
    P = A*P*A'+Q; % cov matrix of the prediction (Pk|k-1)



end

% Estimate + reference display
figure
plot(t,strada.x,t,Xs(1,:)','r');
ylabel('m');
xlabel('t (s)');
title('Estimate and ground truth');
legend('Ground truth','Estimate');

% Errors display with +/- 3 sigma bounds
figure;
plot(t,Xs(1,:)'-strada.x);zoom on;hold on;
plot(t,3*sqrt(Px1),'r');plot(t,-3*sqrt(Px1),'r');ylabel('x error');
xlabel('t (s)');
title('Estimation error with +/- 3 sigma bounds');

disp(['Error mean in x= ', num2str(mean(Xs(1,:)'-strada.x)),...
      '. Error max in x= ', num2str(max(abs(Xs(1,:)'-strada.x)))]);
