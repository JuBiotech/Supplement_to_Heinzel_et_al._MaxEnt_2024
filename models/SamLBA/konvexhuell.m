function [A,center,vol,pax,paxlen] = konvexhuell(dim)
% A := Matrix-Form des in den Ursprung translatierten Ellipsoiden
% center := Zentrum des Ellipsoiden
% vol := Volumen des Ellipsoiden
% pax := Matrix mit Hauptachsen (spaltenweise)
% paxlen := Längen der Hauptachsen

% Exakte Volumenberechnung im Polyeder wird schwierig:
%   http://mathworld.wolfram.com/ConvexPolyhedron.html
%   No method is known for computing the volume of a general convex
%   polyhedron (Ogilvy 1990, p. 173).
%   Ogilvy, C. S. Excursions in Geometry. New York: Dover, 1990.

% % Hyperwürfel-Ecken
% V=zeros(dim,2^dim);
% for k=1:2^dim
% 	V(:,k) = dec2bin(k-1,dim)-'0';
% end
% V = flipud(V);

% Punktwolke
V=2*randn(dim,1000);
% boolean-Vektor (1:=Punkt in Hülle, 0:=Punkt im Polyeder)
inhull = zeros(1,size(V,2));

S = [];
for j=1:2^(dim+3)
    % Zufälliger Richtungsvektor
    d = randvec(dim)';

    % Projektionsebene (senkrecht zur Diagonalen)
    % orthogonale Basis-Vektoren mittels QR
    [Q,R] = qr(d); A = Q(:,2:3);
    
    % Projektion in nd-Ebene mit Basisvektoren A(:,1), A(:,2)
    P = A*inv(A'*A)*A';
    Vp = P*V;
    % Rotation nd-Ebene auf 2d-Ebene
    %Vp2d = pinv(A)*Vp;
    P2d = [qrpkernel(A',[1;0])*[1;ones(dim-2,1)], qrpkernel(A',[0;1])*[1;ones(dim-2,1)]]';
%   disp('probe');
%   P2d * A
    Vp2d = P2d*Vp;

    % konvexe Hülle in der projizierten Punkte berechnen
    hull = zeros(1,size(V,2));
    hull(convhull(Vp2d(1,:),Vp2d(2,:))) = 1;
    inhull = inhull | hull;
    
    S = [ S sum(inhull) ];
end
plot(1:length(S),S);
disp('press enter, dude');
pause;

X = V(:,inhull)';

% Visualisierung für 3d
if dim == 3
    clf;
    hold on;
    plot3(V(1,:),V(2,:),V(3,:),'kx');
    plot3(X(:,1),X(:,2),X(:,3),'ro');
    axis equal;
    hold off;
    disp('press enter, dude');
    pause;
    T=delaunayn(X);
    tetramesh(T,X);
end
[A,center] = MinVolEllipse(X',.01)
[pax,ew] = eig(A);
paxlen = sqrt(1./diag(ew));
vol = ellipsoid_volume(paxlen);

function v = randvec(d)
% Würfelt gleichverteilte Zufallspunkte auf der Oberfläche einer
% Einheitskugel
v = randn(1,d);
r = sqrt(sum(v.*v));
v = v./r;

function [K,p]=qrpkernel(A,b)
% Kernel-Matrix mittels QR(P)-Faktorisierung.
% Zurückgegeben wird die Kernelmatrix K und der Vektor p, der
% die Indizes der freien Variablen enthält.
[Q,R,P] = qr(A);
[m,n] = size(A);
r = rank(A);
if nargin==1
	b = zeros(m,1);
end
Rd = R(1:r,1:r);
invRd = inv(Rd);
Rf = R(1:r,r+1:n);
Q1 = Q(1:m,1:r); % = ran(A)
%Q2 = Q(1:m,r+1:m); % = null(A')
K = P*[ invRd*Q1'*b, -invRd*Rf; zeros(n-r,1), eye(n-r) ];
% Indizes freier Variablen aus P herauskopieren
p = zeros(1,n-r);
for j=r+1:n
	for i=1:n
		if P(i,j)==1
			p(j-r) = i;
			break;
		end
	end
end

function V=ellipsoid_volume(a)
% siehe:
% http://www.ebyte.it/library/docs/math05a/nDimEllipsoidVolumes05.html
nu = @(x)pi^x/gamma(x+1);
n = length(a);
V = nu(n/2) * prod(a);
