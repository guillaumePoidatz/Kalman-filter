% Ph. Bonnifait
clear
clc
load simulated_data.mat; %load variables in the workspace

nepoch=length(t);

%size of the state
n=2; %for instance X=(x,y,vx,vy) - to be changed

% Variables used to store the results
Xs=zeros(n,nepoch); 
%If we want to store some terms of the covariance matrix, wr allocate
%before some memory
Px1=zeros(1,nepoch); %variance
Px1x2=zeros(1,nepoch); %covariance

% Initial state of the filter (to be changed)
X=zeros(n,1);
P=100*eye(n,n);
C = [1,0;0,1];
Qv = 1; % standard deviation of the noise on the model

for i=1:nepoch
    % first we compute the matrices of the model
    Q = [1/4*Qv,0;0,Qv].^2; % random noise on the model
    A = [1,1;0,1]; % no command
    R = [gps.sx,0;0,tachy.sv].^2; % random noise on the measurements
    Y = [gps.x(i);tachy.v(i)];
    % measurement update (estimation)
    % error on the estimation
    epsilon(1,:) = Y(1,:)-C(1,:)*X;
    epsilon(2,:) = Y(2,:)-C(2,:)*X;
    K = P*C*1/(C*P*C'+R); % Kalman gain (Pxy)*Py-1
    %K(2,:) = P(2,:)*C*1/(C*P*C'+R(2,:));
    X(1,:) = X(1,:)+K(1,:)*epsilon; % update of X (Xk|k)
    X(2,:) = X(2,:)+K(2,:)*epsilon;
    P(1,:) = ([1,0]-K(1,:)*C)*P*(1-K*C)'+K(1,:)*R*K'; % cov matrix of the update (Pk|k)
    P(2,:) = ([0,1]-K(2,:)*C)*P*(1-K*C)'+K(2,:)*R*K'; % cov matrix of the update (Pk|k)
    
    % Storage (corresponds to the output of the filter)
    Xs(:,i)=X;
    Px1(i)=P(1,1);
    Px1x2(i)=P(1,2);
    %...
    
    % time update (prediction)
    X(1,:) = A(1,:)*X; % prediction (Xk|k-1)
    X(2,:) = A(2,:)*X;
    % covariance matrix of the estimate
    P(1,:) = A(1,:)*P*A'+ Q(1,:); % cov matrix of the prediction (Pk|k-1)
    P(2,:) = A(2,:)*P*A'+Q(2,:);



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
