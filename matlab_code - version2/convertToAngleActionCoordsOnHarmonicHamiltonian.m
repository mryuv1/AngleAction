function AAVars = convertToAngleActionCoordsOnHarmonicHamiltonian(q, w) 
% q = [x, y, px, py]; w = [w1, w2]

x = q(:,1);
y = q(:,2);
px = q(:,3);
py = q(:,4);
w1 = w(1);
w2 = w(2);

% with harmonic Hamiltonian we can define:
% q = sqrt(2I/w)cos(theta), p = -sqrt(Iw)sin(theta)
% from that we can extract:
% theta = atan(qmw/p), J = E/w

Ex = 0.5*(px.^2) + 0.5*((w1*x).^2);
Ey = 0.5*(py.^2) + 0.5*((w2*y).^2);

Jx = Ex./w1;
Jy = Ey./w2;


Thetax = atan2(-px, x*w1);
Thetay = atan2(-py, y*w2);

Thetax(abs(px)<1e-12)=0;
Thetay(abs(py)<1e-12)=0;

AAVars = struct();
AAVars.Jx = Jx;
AAVars.Jy = Jy;
AAVars.Thetax = Thetax;
AAVars.Thetay = Thetay;