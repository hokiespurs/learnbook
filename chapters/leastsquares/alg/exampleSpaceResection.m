function exampleSpaceResection
%% Test using least squares for a single photo resection
%Assumes No Error in Camera Interior Orientation

%% User Inputs
RANDOMSEED = 1;
% camera position
CAM.TRUEXYZ = [0 0 100];
CAM.TRUERPYDEG = [10 180 20];
% least squares initialization guess
LSQINITIALXYZ = [0 0 100];
LSQINITIALRPY = [0 180 0];
% camera parameters
CAM.XPIX = 640;
CAM.YPIX = 480;
CAM.HFOV = 90;
CAM.CX = 320;
CAM.CY = 240;
% add GCPs
NGCP = 6;
GCPX = [-50 50];
GCPY = [-50 50];
GCPZ = [-20 20];
% Add Noise to simulated data?
PIXNOISE = 1;
COORDNOISE = 0.1;
% Estimate of Noise into Total Least Squares
PIXNOISEEST = 1;
COORDNOISEEST = 0.1;

%% Generate Simulated Data for Control Points and Pixel Coordinates
rng(RANDOMSEED);
% Make K Matrix
f = CAM.XPIX/2/tand(CAM.HFOV/2);
K = [f 0 CAM.CX;0 f CAM.CY;0 0 1];
CAM.TRUERPYRAD = CAM.TRUERPYDEG*pi/180;
% Make R matrix
R = rpy2dcm(CAM.TRUERPYRAD(1),CAM.TRUERPYRAD(2),CAM.TRUERPYRAD(3));
% Calculate Sample GCPs
truegcpx = GCPX(1)+rand(NGCP,1)*(GCPX(2)-GCPX(1));
truegcpy = GCPY(1)+rand(NGCP,1)*(GCPY(2)-GCPY(1));
truegcpz = GCPZ(1)+rand(NGCP,1)*(GCPZ(2)-GCPZ(1));
% project into uv space
uvq = K*R*[truegcpx' - CAM.TRUEXYZ(1);...
           truegcpy' - CAM.TRUEXYZ(2);...
           truegcpz' - CAM.TRUEXYZ(3)];
trueu = uvq(1,:)'./uvq(3,:)';
truev = uvq(2,:)'./uvq(3,:)';
%% Add Noise to Simulated Data
noisygcpx = truegcpx + randn(numel(truegcpx),1)*COORDNOISE;
noisygcpy = truegcpy + randn(numel(truegcpx),1)*COORDNOISE;
noisygcpz = truegcpz + randn(numel(truegcpx),1)*COORDNOISE;
noisyu = trueu + randn(numel(truegcpx),1)*PIXNOISE;
noisyv = truev + randn(numel(truegcpx),1)*PIXNOISE;
%% Visualize Data
plotCamAndGCPs(K,R,CAM.TRUEXYZ',[truegcpx,truegcpy,truegcpz]',CAM.XPIX,CAM.YPIX);

%% Generate Covariance Matrix for Observation Values
covObs = kron(eye(NGCP),...
    diag([ones(1,2)*PIXNOISEEST,ones(1,3)*COORDNOISEEST]));

%% Anonymous Functions 
Jfun = @(Xi) Jfunction(Xi(1),Xi(2),Xi(3),Xi(4),Xi(5),Xi(6),...
	noisyu,noisyv,noisygcpx,noisygcpy,noisygcpz,f,CAM.CX,CAM.CY,...
	 logical([1 1 1 1 1 1 ]));
Bfun = @(Xi) Bfunction(Xi(1),Xi(2),Xi(3),Xi(4),Xi(5),Xi(6),...
	noisyu,noisyv,noisygcpx,noisygcpy,noisygcpz,f,CAM.CX,CAM.CY);
Kfun = @(Xi) Kfunction(Xi(1),Xi(2),Xi(3),Xi(4),Xi(5),Xi(6),...
	noisyu,noisyv,noisygcpx,noisygcpy,noisygcpz,f,CAM.CX,CAM.CY);

%% Set Initial Least Squares Guess Xo
Xo = [LSQINITIALRPY*pi/180 LSQINITIALXYZ]';

%% Perform Total Least Squares
fprintf('Performing Total Least Squares\n');
[X_tls,Sx_tls,lsrinfo_tls] = lsrtls(Jfun,Kfun,Bfun,Xo,covObs);
covXYZ = enforceCovariance(Sx_tls(4:6,4:6));
plotCovarianceSurf(X_tls(4),X_tls(5),X_tls(6),covXYZ,0.95,lsrinfo_tls.dof);
%% Perform Nonlinear Least Squares
fprintf('Performing Unqeighted Nonlinear Least Squares\n');
[X_nlin,Sx_nlin,lsrinfo_nlin] = lsrnlin(Jfun,Kfun,Xo);

%% Perform TLS with Z Elevation Locked at True Z Position
fprintf('Performing Total Least Squares with Z Locked\n');

% Set New Anonymous Functions 
Jfun = @(Xi) Jfunction(Xi(1),Xi(2),Xi(3),Xi(4),Xi(5),CAM.TRUEXYZ(3),...
	noisyu,noisyv,noisygcpx,noisygcpy,noisygcpz,f,CAM.CX,CAM.CY,...
	 logical([1 1 1 1 1 0]));
Bfun = @(Xi) Bfunction(Xi(1),Xi(2),Xi(3),Xi(4),Xi(5),CAM.TRUEXYZ(3),...
	noisyu,noisyv,noisygcpx,noisygcpy,noisygcpz,f,CAM.CX,CAM.CY);
Kfun = @(Xi) Kfunction(Xi(1),Xi(2),Xi(3),Xi(4),Xi(5),CAM.TRUEXYZ(3),...
	noisyu,noisyv,noisygcpx,noisygcpy,noisygcpz,f,CAM.CX,CAM.CY);
 
 % need a new Xo, since not solving for Z
 Xo = [LSQINITIALRPY*pi/180 LSQINITIALXYZ(1:2)]';

 [X_tlsLockZ,Sx_tlsLockZ,lsrinfo_tlsLockZ] = lsrtls(Jfun,Kfun,Bfun,Xo,covObs);
 end
function K = Kfunction(roll,pitch,yaw,Xc,Yc,Zc,u,v,Xw,Yw,Zw,f,cx,cy)
%% K EQUATIONS 
K = nan(numel(u)*2,1);
K(1:2:end,1) = u + ((cx.*(cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)) - f.*cos(pitch).*sin(yaw)).*(Yc - Yw) + (Zc - Zw).*(f.*sin(pitch) - cx.*cos(pitch).*cos(roll)) - (cx.*(sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch)) + f.*cos(pitch).*cos(yaw)).*(Xc - Xw))./((Xc - Xw).*(sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch)) - (Yc - Yw).*(cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)) + cos(pitch).*cos(roll).*(Zc - Zw));
K(2:2:end,1) = v - ((cy.*cos(pitch).*cos(roll) + f.*cos(pitch).*sin(roll)).*(Zc - Zw) + (Xc - Xw).*(cy.*(sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch)) - f.*(cos(roll).*sin(yaw) - cos(yaw).*sin(pitch).*sin(roll))) - (Yc - Yw).*(cy.*(cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)) - f.*(cos(roll).*cos(yaw) + sin(pitch).*sin(roll).*sin(yaw))))./((Xc - Xw).*(sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch)) - (Yc - Yw).*(cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)) + cos(pitch).*cos(roll).*(Zc - Zw));
end

function J = Jfunction(roll,pitch,yaw,Xc,Yc,Zc,u,v,Xw,Yw,Zw,f,cx,cy,flag)
%% J EQUATIONS 
J = nan(numel(u)*2,6);
J(1:2:end,1) = - (cx.*(Yc - Yw).*(cos(roll).*cos(yaw) + sin(pitch).*sin(roll).*sin(yaw)) - cx.*(Xc - Xw).*(cos(roll).*sin(yaw) - cos(yaw).*sin(pitch).*sin(roll)) + cx.*cos(pitch).*sin(roll).*(Zc - Zw))./((Xc - Xw).*(sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch)) - (Yc - Yw).*(cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)) + cos(pitch).*cos(roll).*(Zc - Zw)) - (((cx.*(cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)) - f.*cos(pitch).*sin(yaw)).*(Yc - Yw) + (Zc - Zw).*(f.*sin(pitch) - cx.*cos(pitch).*cos(roll)) - (cx.*(sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch)) + f.*cos(pitch).*cos(yaw)).*(Xc - Xw)).*((Yc - Yw).*(cos(roll).*cos(yaw) + sin(pitch).*sin(roll).*sin(yaw)) - (Xc - Xw).*(cos(roll).*sin(yaw) - cos(yaw).*sin(pitch).*sin(roll)) + cos(pitch).*sin(roll).*(Zc - Zw)))./((Xc - Xw).*(sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch)) - (Yc - Yw).*(cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)) + cos(pitch).*cos(roll).*(Zc - Zw)).^2;
J(1:2:end,2) = (((cx.*(cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)) - f.*cos(pitch).*sin(yaw)).*(Yc - Yw) + (Zc - Zw).*(f.*sin(pitch) - cx.*cos(pitch).*cos(roll)) - (cx.*(sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch)) + f.*cos(pitch).*cos(yaw)).*(Xc - Xw)).*(cos(pitch).*cos(roll).*cos(yaw).*(Xc - Xw) - cos(roll).*sin(pitch).*(Zc - Zw) + cos(pitch).*cos(roll).*sin(yaw).*(Yc - Yw)))./((Xc - Xw).*(sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch)) - (Yc - Yw).*(cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)) + cos(pitch).*cos(roll).*(Zc - Zw)).^2 - ((f.*sin(pitch).*sin(yaw) - cx.*cos(pitch).*cos(roll).*sin(yaw)).*(Yc - Yw) + (Zc - Zw).*(f.*cos(pitch) + cx.*cos(roll).*sin(pitch)) + (f.*cos(yaw).*sin(pitch) - cx.*cos(pitch).*cos(roll).*cos(yaw)).*(Xc - Xw))./((Xc - Xw).*(sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch)) - (Yc - Yw).*(cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)) + cos(pitch).*cos(roll).*(Zc - Zw));
J(1:2:end,3) = ((cx.*(sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch)) + f.*cos(pitch).*cos(yaw)).*(Yc - Yw) + (cx.*(cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)) - f.*cos(pitch).*sin(yaw)).*(Xc - Xw))./((Xc - Xw).*(sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch)) - (Yc - Yw).*(cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)) + cos(pitch).*cos(roll).*(Zc - Zw)) + (((Xc - Xw).*(cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)) + (Yc - Yw).*(sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch))).*((cx.*(cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)) - f.*cos(pitch).*sin(yaw)).*(Yc - Yw) + (Zc - Zw).*(f.*sin(pitch) - cx.*cos(pitch).*cos(roll)) - (cx.*(sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch)) + f.*cos(pitch).*cos(yaw)).*(Xc - Xw)))./((Xc - Xw).*(sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch)) - (Yc - Yw).*(cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)) + cos(pitch).*cos(roll).*(Zc - Zw)).^2;
J(1:2:end,4) = (cx.*(sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch)) + f.*cos(pitch).*cos(yaw))./((Xc - Xw).*(sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch)) - (Yc - Yw).*(cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)) + cos(pitch).*cos(roll).*(Zc - Zw)) + ((sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch)).*((cx.*(cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)) - f.*cos(pitch).*sin(yaw)).*(Yc - Yw) + (Zc - Zw).*(f.*sin(pitch) - cx.*cos(pitch).*cos(roll)) - (cx.*(sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch)) + f.*cos(pitch).*cos(yaw)).*(Xc - Xw)))./((Xc - Xw).*(sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch)) - (Yc - Yw).*(cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)) + cos(pitch).*cos(roll).*(Zc - Zw)).^2;
J(1:2:end,5) = - (cx.*(cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)) - f.*cos(pitch).*sin(yaw))./((Xc - Xw).*(sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch)) - (Yc - Yw).*(cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)) + cos(pitch).*cos(roll).*(Zc - Zw)) - ((cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)).*((cx.*(cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)) - f.*cos(pitch).*sin(yaw)).*(Yc - Yw) + (Zc - Zw).*(f.*sin(pitch) - cx.*cos(pitch).*cos(roll)) - (cx.*(sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch)) + f.*cos(pitch).*cos(yaw)).*(Xc - Xw)))./((Xc - Xw).*(sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch)) - (Yc - Yw).*(cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)) + cos(pitch).*cos(roll).*(Zc - Zw)).^2;
J(1:2:end,6) = (cos(pitch).*cos(roll).*((cx.*(cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)) - f.*cos(pitch).*sin(yaw)).*(Yc - Yw) + (Zc - Zw).*(f.*sin(pitch) - cx.*cos(pitch).*cos(roll)) - (cx.*(sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch)) + f.*cos(pitch).*cos(yaw)).*(Xc - Xw)))./((Xc - Xw).*(sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch)) - (Yc - Yw).*(cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)) + cos(pitch).*cos(roll).*(Zc - Zw)).^2 - (f.*sin(pitch) - cx.*cos(pitch).*cos(roll))./((Xc - Xw).*(sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch)) - (Yc - Yw).*(cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)) + cos(pitch).*cos(roll).*(Zc - Zw));
J(2:2:end,1) = ((f.*cos(pitch).*cos(roll) - cy.*cos(pitch).*sin(roll)).*(Zc - Zw) + (Xc - Xw).*(cy.*(cos(roll).*sin(yaw) - cos(yaw).*sin(pitch).*sin(roll)) + f.*(sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch))) - (Yc - Yw).*(cy.*(cos(roll).*cos(yaw) + sin(pitch).*sin(roll).*sin(yaw)) + f.*(cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw))))./((Xc - Xw).*(sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch)) - (Yc - Yw).*(cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)) + cos(pitch).*cos(roll).*(Zc - Zw)) + (((cy.*cos(pitch).*cos(roll) + f.*cos(pitch).*sin(roll)).*(Zc - Zw) + (Xc - Xw).*(cy.*(sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch)) - f.*(cos(roll).*sin(yaw) - cos(yaw).*sin(pitch).*sin(roll))) - (Yc - Yw).*(cy.*(cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)) - f.*(cos(roll).*cos(yaw) + sin(pitch).*sin(roll).*sin(yaw)))).*((Yc - Yw).*(cos(roll).*cos(yaw) + sin(pitch).*sin(roll).*sin(yaw)) - (Xc - Xw).*(cos(roll).*sin(yaw) - cos(yaw).*sin(pitch).*sin(roll)) + cos(pitch).*sin(roll).*(Zc - Zw)))./((Xc - Xw).*(sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch)) - (Yc - Yw).*(cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)) + cos(pitch).*cos(roll).*(Zc - Zw)).^2;
J(2:2:end,2) = (f.*(Yc.*cos(yaw) - Yw.*cos(yaw) - Xc.*sin(yaw) + Xw.*sin(yaw)).*(Zc.*sin(pitch) - Zw.*sin(pitch) - Xc.*cos(pitch).*cos(yaw) + Xw.*cos(pitch).*cos(yaw) - Yc.*cos(pitch).*sin(yaw) + Yw.*cos(pitch).*sin(yaw)))./(Zc.*cos(pitch).*cos(roll) - Zw.*cos(pitch).*cos(roll) - Yc.*cos(yaw).*sin(roll) + Yw.*cos(yaw).*sin(roll) + Xc.*sin(roll).*sin(yaw) - Xw.*sin(roll).*sin(yaw) + Xc.*cos(roll).*cos(yaw).*sin(pitch) - Xw.*cos(roll).*cos(yaw).*sin(pitch) + Yc.*cos(roll).*sin(pitch).*sin(yaw) - Yw.*cos(roll).*sin(pitch).*sin(yaw)).^2;
J(2:2:end,3) = ((Xc - Xw).*(cy.*(cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)) - f.*(cos(roll).*cos(yaw) + sin(pitch).*sin(roll).*sin(yaw))) + (Yc - Yw).*(cy.*(sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch)) - f.*(cos(roll).*sin(yaw) - cos(yaw).*sin(pitch).*sin(roll))))./((Xc - Xw).*(sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch)) - (Yc - Yw).*(cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)) + cos(pitch).*cos(roll).*(Zc - Zw)) - (((Xc - Xw).*(cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)) + (Yc - Yw).*(sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch))).*((cy.*cos(pitch).*cos(roll) + f.*cos(pitch).*sin(roll)).*(Zc - Zw) + (Xc - Xw).*(cy.*(sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch)) - f.*(cos(roll).*sin(yaw) - cos(yaw).*sin(pitch).*sin(roll))) - (Yc - Yw).*(cy.*(cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)) - f.*(cos(roll).*cos(yaw) + sin(pitch).*sin(roll).*sin(yaw)))))./((Xc - Xw).*(sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch)) - (Yc - Yw).*(cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)) + cos(pitch).*cos(roll).*(Zc - Zw)).^2;
J(2:2:end,4) = -(f.*(Yc.*sin(pitch) - Yw.*sin(pitch) + Zc.*cos(pitch).*sin(yaw) - Zw.*cos(pitch).*sin(yaw)))./(Zc.*cos(pitch).*cos(roll) - Zw.*cos(pitch).*cos(roll) - Yc.*cos(yaw).*sin(roll) + Yw.*cos(yaw).*sin(roll) + Xc.*sin(roll).*sin(yaw) - Xw.*sin(roll).*sin(yaw) + Xc.*cos(roll).*cos(yaw).*sin(pitch) - Xw.*cos(roll).*cos(yaw).*sin(pitch) + Yc.*cos(roll).*sin(pitch).*sin(yaw) - Yw.*cos(roll).*sin(pitch).*sin(yaw)).^2;
J(2:2:end,5) = (f.*(Xc.*sin(pitch) - Xw.*sin(pitch) + Zc.*cos(pitch).*cos(yaw) - Zw.*cos(pitch).*cos(yaw)))./(Zc.*cos(pitch).*cos(roll) - Zw.*cos(pitch).*cos(roll) - Yc.*cos(yaw).*sin(roll) + Yw.*cos(yaw).*sin(roll) + Xc.*sin(roll).*sin(yaw) - Xw.*sin(roll).*sin(yaw) + Xc.*cos(roll).*cos(yaw).*sin(pitch) - Xw.*cos(roll).*cos(yaw).*sin(pitch) + Yc.*cos(roll).*sin(pitch).*sin(yaw) - Yw.*cos(roll).*sin(pitch).*sin(yaw)).^2;
J(2:2:end,6) = -(f.*cos(pitch).*(Yc.*cos(yaw) - Yw.*cos(yaw) - Xc.*sin(yaw) + Xw.*sin(yaw)))./(Zc.*cos(pitch).*cos(roll) - Zw.*cos(pitch).*cos(roll) - Yc.*cos(yaw).*sin(roll) + Yw.*cos(yaw).*sin(roll) + Xc.*sin(roll).*sin(yaw) - Xw.*sin(roll).*sin(yaw) + Xc.*cos(roll).*cos(yaw).*sin(pitch) - Xw.*cos(roll).*cos(yaw).*sin(pitch) + Yc.*cos(roll).*sin(pitch).*sin(yaw) - Yw.*cos(roll).*sin(pitch).*sin(yaw)).^2;
J(:,~flag)=[]; % handle optional unknowns
end

function B = Bfunction(roll,pitch,yaw,Xc,Yc,Zc,u,v,Xw,Yw,Zw,f,cx,cy)
%% B EQUATIONS 
B = nan(numel(u)*2,5);
B(1:2:end,1) = -1;
B(1:2:end,2) = 0;
B(1:2:end,3) = - (cx.*(sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch)) + f.*cos(pitch).*cos(yaw))./((Xc - Xw).*(sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch)) - (Yc - Yw).*(cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)) + cos(pitch).*cos(roll).*(Zc - Zw)) - ((sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch)).*((cx.*(cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)) - f.*cos(pitch).*sin(yaw)).*(Yc - Yw) + (Zc - Zw).*(f.*sin(pitch) - cx.*cos(pitch).*cos(roll)) - (cx.*(sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch)) + f.*cos(pitch).*cos(yaw)).*(Xc - Xw)))./((Xc - Xw).*(sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch)) - (Yc - Yw).*(cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)) + cos(pitch).*cos(roll).*(Zc - Zw)).^2;
B(1:2:end,4) = (cx.*(cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)) - f.*cos(pitch).*sin(yaw))./((Xc - Xw).*(sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch)) - (Yc - Yw).*(cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)) + cos(pitch).*cos(roll).*(Zc - Zw)) + ((cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)).*((cx.*(cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)) - f.*cos(pitch).*sin(yaw)).*(Yc - Yw) + (Zc - Zw).*(f.*sin(pitch) - cx.*cos(pitch).*cos(roll)) - (cx.*(sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch)) + f.*cos(pitch).*cos(yaw)).*(Xc - Xw)))./((Xc - Xw).*(sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch)) - (Yc - Yw).*(cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)) + cos(pitch).*cos(roll).*(Zc - Zw)).^2;
B(1:2:end,5) = (f.*sin(pitch) - cx.*cos(pitch).*cos(roll))./((Xc - Xw).*(sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch)) - (Yc - Yw).*(cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)) + cos(pitch).*cos(roll).*(Zc - Zw)) - (cos(pitch).*cos(roll).*((cx.*(cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)) - f.*cos(pitch).*sin(yaw)).*(Yc - Yw) + (Zc - Zw).*(f.*sin(pitch) - cx.*cos(pitch).*cos(roll)) - (cx.*(sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch)) + f.*cos(pitch).*cos(yaw)).*(Xc - Xw)))./((Xc - Xw).*(sin(roll).*sin(yaw) + cos(roll).*cos(yaw).*sin(pitch)) - (Yc - Yw).*(cos(yaw).*sin(roll) - cos(roll).*sin(pitch).*sin(yaw)) + cos(pitch).*cos(roll).*(Zc - Zw)).^2;
B(2:2:end,1) = 0;
B(2:2:end,2) = -1;
B(2:2:end,3) = (f.*(Yc.*sin(pitch) - Yw.*sin(pitch) + Zc.*cos(pitch).*sin(yaw) - Zw.*cos(pitch).*sin(yaw)))./(Zc.*cos(pitch).*cos(roll) - Zw.*cos(pitch).*cos(roll) - Yc.*cos(yaw).*sin(roll) + Yw.*cos(yaw).*sin(roll) + Xc.*sin(roll).*sin(yaw) - Xw.*sin(roll).*sin(yaw) + Xc.*cos(roll).*cos(yaw).*sin(pitch) - Xw.*cos(roll).*cos(yaw).*sin(pitch) + Yc.*cos(roll).*sin(pitch).*sin(yaw) - Yw.*cos(roll).*sin(pitch).*sin(yaw)).^2;
B(2:2:end,4) = -(f.*(Xc.*sin(pitch) - Xw.*sin(pitch) + Zc.*cos(pitch).*cos(yaw) - Zw.*cos(pitch).*cos(yaw)))./(Zc.*cos(pitch).*cos(roll) - Zw.*cos(pitch).*cos(roll) - Yc.*cos(yaw).*sin(roll) + Yw.*cos(yaw).*sin(roll) + Xc.*sin(roll).*sin(yaw) - Xw.*sin(roll).*sin(yaw) + Xc.*cos(roll).*cos(yaw).*sin(pitch) - Xw.*cos(roll).*cos(yaw).*sin(pitch) + Yc.*cos(roll).*sin(pitch).*sin(yaw) - Yw.*cos(roll).*sin(pitch).*sin(yaw)).^2;
B(2:2:end,5) = (f.*cos(pitch).*(Yc.*cos(yaw) - Yw.*cos(yaw) - Xc.*sin(yaw) + Xw.*sin(yaw)))./(Zc.*cos(pitch).*cos(roll) - Zw.*cos(pitch).*cos(roll) - Yc.*cos(yaw).*sin(roll) + Yw.*cos(yaw).*sin(roll) + Xc.*sin(roll).*sin(yaw) - Xw.*sin(roll).*sin(yaw) + Xc.*cos(roll).*cos(yaw).*sin(pitch) - Xw.*cos(roll).*cos(yaw).*sin(pitch) + Yc.*cos(roll).*sin(pitch).*sin(yaw) - Yw.*cos(roll).*sin(pitch).*sin(yaw)).^2;
B = bumphdiag(B,2); %this pads the matrix with 0s
end

function R = rpy2dcm(roll,pitch,yaw)
% Calculate Direction Cosine Matrix Using Tait-Bryan (Z-Y-X convention)
Rx = [1 0 0; 0 cos(roll) sin(roll); 0 -sin(roll) cos(roll)];
Ry = [cos(pitch) 0 -sin(pitch); 0 1 0;sin(pitch) 0 cos(pitch)];
Rz = [cos(yaw) sin(yaw) 0; -sin(yaw) cos(yaw) 0; 0 0 1];
R = Rx*Ry*Rz;
end

function plotCamAndGCPs(K,R,camT,gcpxyz,xpix,ypix)
% this code is kinda old and probably uses some sketchy backwards math...
% but it makes a cool figure!

% calculate pixel coordinates
iRT = inv([inv(R) camT;0 0 0 1]);
iT = iRT(1:3,4);

uv1 = K * [R iT] * [gcpxyz; ones(size(gcpxyz(1,:)))];
u = uv1(1,:)./uv1(3,:);
v = uv1(2,:)./uv1(3,:);

maxdist = max(uv1(3,:));

FSCALE = K(1,1)/xpix;
CAMSCALE = maxdist/4/FSCALE;

%cam fov points
xc = CAMSCALE * [0 -1 -1 1 1 -1 1 0 1 -1 0]*0.5;
yc = CAMSCALE * ypix/xpix * [0 1 -1 -1 1 1 1 0 -1 -1 0]*0.5;
zc = CAMSCALE * FSCALE * [0 1 1 1 1 1 1 0 1 1 0];
xyzc = inv(R) * [xc; yc; zc];

% cam u,v points
uplanecoord = CAMSCALE * (u - K(1,3))/ xpix ;
vplanecoord = CAMSCALE * ypix/xpix * (v- K(2,3)) / ypix ;
zplanecoord = CAMSCALE * FSCALE * ones(size(u));

planecoordxyz = inv(R)*[uplanecoord; vplanecoord;zplanecoord];

% cam rectangle plane
ind = 2:5;
xp = xyzc(1,ind);
yp = xyzc(2,ind);
zp = xyzc(3,ind);

figure(1);clf
%plot camera fov
plot3(camT(1) + xyzc(1,:),camT(2) + xyzc(2,:),camT(3) + xyzc(3,:),'k')
hold on

%plot camera fov plane
patch(camT(1) + xp,camT(2) + yp,camT(3) + zp,'k','FaceAlpha',0.2)

%plot camera to gcp lines
for i=1:numel(gcpxyz(1,:))
   % gcp lines
   plot3([camT(1) gcpxyz(1,i)],[camT(2) gcpxyz(2,i)],[camT(3) gcpxyz(3,i)],'r.-','markersize',20) 
end
% uv coords
plot3(camT(1) + planecoordxyz(1,:),camT(2) + planecoordxyz(2,:),camT(3) + planecoordxyz(3,:),'g.','markersize',20)
indbad = u > xpix | u<0 | v>ypix |v<0 | uv1(3,:)<0;
plot3(camT(1) + planecoordxyz(1,indbad),camT(2) + planecoordxyz(2,indbad),camT(3) + planecoordxyz(3,indbad),'m.','markersize',20)


% plot camera origin coordinate
plot3(camT(1),camT(2),camT(3),'b.','markersize',20);

xlabel('X Coordinate','fontsize',20);
ylabel('Y Coordinate','fontsize',20);
zlabel('Z Coordinate','fontsize',20)
axis equal
end

function symbolicKJBequations
%% Generate J,B,K functions for space resection equations using symbolic toolbox
syms u v Xw Yw Zw roll pitch yaw Xc Yc Zc f cx cy
% Calculate Direction Cosine Matrix Using Tait-Bryan (Z-Y-X convention)
Rx = [1 0 0; 0 cos(roll) sin(roll); 0 -sin(roll) cos(roll)];
Ry = [cos(pitch) 0 -sin(pitch); 0 1 0;sin(pitch) 0 cos(pitch)];
Rz = [cos(yaw) sin(yaw) 0; -sin(yaw) cos(yaw) 0; 0 0 1];
R = Rx*Ry*Rz;

uvq =[f 0 cx;0 f cy;0 0 1]*R*[Xw-Xc;Yw-Yc;Zw-Zc];

F = [u - uvq(1)/uvq(3);v - uvq(2)/uvq(3)];

betavars = [roll pitch yaw Xc Yc Zc];
invars = [u v Xw Yw Zw];
[J,B,K]=calcKJBequations(F,invars,betavars);

end

function hx = bumphdiag(x,n)
hx = [];
for i=n:n:size(x,1)
    hx = blkdiag(hx,x(i-n+1:i,:));
end

end