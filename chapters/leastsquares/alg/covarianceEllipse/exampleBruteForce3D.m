%%
Cxx = 1;Cyy = 2;Czz = 4;
Cxy = 1;Cxz = .5;Cyz = 2;
%MAG
Cxx = 1.394500645924E-005;
Cyy = 2.371056244075E-005;
Czz = 9.695190117622E-005;
Cxy = 9.368919368385E-006;
Cxz = -2.005312258009E-005;
Cyz = -8.209277685364E-006;

C = [Cxx Cxy Cxz;
      Cxy Cyy Cyz;
      Cxz Cyz Czz];

%
[V,D]=eig(C);
[az,el,r]=cart2sph(V(1,:),V(2,:),V(3,:));
eValue = sqrt(diag(D));

ielevs = linspace(-pi/2,pi/2,50);
iazs = linspace(-pi,pi,50);

[u,v]=meshgrid(ielevs,iazs);
x = eValue(1).*cos(u).*cos(v);
y = eValue(2).*cos(u).*sin(v);
z = eValue(3).*sin(u);
r = sqrt((x).^2+(y).^2+(z).^2);

xyz = [x(:) y(:) z(:)];
xyz2 = V * xyz';

x = reshape(xyz2(1,:),size(x));
y = reshape(xyz2(2,:),size(y));
z = reshape(xyz2(3,:),size(z));

figure(2);clf
maxval = max(eValue);
% Plot 3D
axes('Position',[0.05 0.35 0.9 0.6])
plotRect([0 0 0],sqrt([Cxx Cyy Czz]),'k.-','linewidth',5); axis equal
hold on
axis equal
surf(x,y,z,r)
alpha(0.3)
hold on
for i=1:3
plot3([0 V(1,i)].*eValue(i),[0 V(2,i)].*eValue(i),[0 V(3,i)].*eValue(i),'linewidth',3)
end
xlabel('X');
ylabel('Y');
zlabel('Z');
axis([-maxval maxval -maxval maxval -maxval maxval])

% Plot X,Z
axes('Position',[0.05 0.05 0.25 0.25])
plotRect([0 0 0],sqrt([Cxx Cyy Czz]),'k.-','linewidth',5); axis equal
hold on
axis equal
surf(x,y,z,r)
alpha(0.3)
hold on
for i=1:3
plot3([0 V(1,i)].*eValue(i),[0 V(2,i)].*eValue(i),[0 V(3,i)].*eValue(i),'linewidth',3)
end
view(0,0);
xlabel('X');
ylabel('Y');
zlabel('Z');
axis([-maxval maxval -maxval maxval -maxval maxval])

% PLOT Y,Z
axes('Position',[0.35 0.05 0.25 0.25])
plotRect([0 0 0],sqrt([Cxx Cyy Czz]),'k.-','linewidth',5); axis equal
hold on
h = surf(x,y,z,r)
alpha(0.3)
hold on
for i=1:3
plot3([0 V(1,i)].*eValue(i),[0 V(2,i)].*eValue(i),[0 V(3,i)].*eValue(i),'linewidth',3)
end
view(90,0)
xlabel('X');
ylabel('Y');
zlabel('Z');
axis([-maxval maxval -maxval maxval -maxval maxval])

% PLOT X,Y
axes('Position',[0.65 0.05 0.25 0.25])
plotRect([0 0 0],sqrt([Cxx Cyy Czz]),'k.-','linewidth',5); axis equal
hold on
surf(x,y,z,r)
alpha(0.3)
hold on
for i=1:3
plot3([0 V(1,i)].*eValue(i),[0 V(2,i)].*eValue(i),[0 V(3,i)].*eValue(i),'linewidth',3)
end
view(0,90)
xlabel('X');
ylabel('Y');
zlabel('Z');
axis([-maxval maxval -maxval maxval -maxval maxval])
