%PROGRAMA QUE OBTIENE DERIVADAS ESPACIALES DE UNA MALLA DE DATOS
%SE NECESITA MALLA EN DOMINIO TIEMPO g(xx,yy) Y MUESTREO EN AMBAS
%DIRECCIONES DX y DY
function[dxg dyg dzg]= DERIVADAS_1(g,dx,dy)
%da las derivadas y recibe malla de datos, dx y dy
[Ny Nx]= size(g);
%CREAR OPERADOR DIFERENCIAL RESPECTO X
DX=[1 0 -1]*(1/(2*dx));
%MAPEO SEÑAL SOBRE OTRA = CONVOLUCIÓN
%OPERADOR DERIVADA RESPECTO Y
DY = [1;0;-1]*(1/(2*dy));
%OPERADORES SE VAN A CONVOLUCIONAR PARA OBTENER DERIVADA:
dxg = conv2(g,DX,'same');
dyg = conv2(g,DY,'same');
%operador derivada se opera en malla de datos, se tendrán problemas en
%frontera
%EDITANDO LAS FRONTERAS
%se copian fronteras en otra columna
dxg(:,1)=dxg(:,2);
dxg(:,Nx)=dxg(:,Nx-1);
dyg(1,:)=dyg(2,:);
dyg(Ny,:)=dyg(Ny-1,:);

%DERIVAR ESPECTRALMENTE EN DOMINIO FOURIER
%CREAR DOMINO FOURIER
FNX = 1/(2*dx);
FNY = 1/(2*dy);
F0X = 1/(Nx*dx);
F0Y = 1/(Ny*dy);
%CREAR VECTORES DIRECCIONALES EN DOMINIO FOURIER
p = [-FNX:F0X:FNX-F0X]*2*pi;
q = [-FNY:F0Y:FNY-F0Y]*2*pi;
%CREANDO EL DOMINIO DE FOURIER
[P Q] = meshgrid(p,q);

%DEFINIR DERIVADA RESPECTO A Z
%APLICANDO DERIVADA EN DOMINIO DE FOURIER DZF(P,Q)=SQRT(P.^2+Q.^2)*G(P,Q)
K = sqrt(P.^2+Q.^2);
%TRANSFORMAR DATOS:
G = fftshift((fft2(g)));
DZG =K.*G;

%REGRESANDO A DOMINIO TIEMPO
dzg = real(ifft2(fftshift(DZG)));