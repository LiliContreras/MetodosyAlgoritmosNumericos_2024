%SE CARGAN LOS DATOS
fx = load('Gx_prisma2.txt');
fy = load('Gy_prisma2.txt');
fz = load('Gz_prisma2.txt');

%generar dominio enteros para realizar Ec. de Euler (x1,x2,...xn y Y)
%numero de observaciones
[Ny Nx]= size(fz);
%vector de posición x
x = [0:Nx-1];
y = [0:Ny-1];
%dominio de enteros
%en función vectores direccionales:
[X Y] = meshgrid(x,y);
%grafica datos obs
figure(1)
mesh(X,Y,fz);
view(2)
title('Anomalía en Dominio de Enteros');
xlabel('Distancia x [km]');
ylabel('Distancia y [km]');
zlabel('Anomalía [mGal]');
c = colorbar;
c.Label.String = '[mGal]';

%SE DEFINEN PARÁMETROS DE EXPLORACIÓN:
%Tolerancia:
Tol = 50;
%Tol = 30;
%Tol = 60;
%Tol = 70;
%Indice Estructural
nu = -1; %PARA UN DIQUE
%Tamaño de Ventana:
v = 20;
%v = 10;
%v = 30;
%v = 40;

%PEDIR MEMORIA PARA ALMACENAR VENTANAS MUESTREADAS:
%ventana asociada a datos observados
Vfx = zeros(v);
Vfy = zeros(v);
Vfz = zeros(v);
%matriz para indice estructural
IE = ones(v)*nu;

%VENTANA DE DERIVADAS DE X y Y
%FX:
Vdxfx = zeros(v);
Vdyfx = zeros(v);
Vdzfx = zeros(v);
%FY:
Vdxfy = zeros(v);
Vdyfy = zeros(v);
Vdzfy = zeros(v);
%FZ:
Vdxfz = zeros(v);
Vdyfz = zeros(v);
Vdzfz = zeros(v);

%VENTANA DOMINIO X y Y
VX = zeros(v);
VY = zeros(v);

%GENERAR MALLAS NECESARIAS PARA CONSTRUIR EC. EULER
dx = 10/(Nx-1); %10 km observación
dy = 10/(Ny-1);
[dxfx dyfx dzfx] = DERIVADAS_1(fx,dx,dy);
[dxfy dyfy dzfy] = DERIVADAS_1(fy,dx,dy);
[dxfz dyfz dzfz] = DERIVADAS_1(fz,dx,dy);

%DOBLE CICLO FOR PARA EXPLORAR DOMINIO CON VENTANAS CUADRADAS DE TAMAÑO V
kx=0; %contador que almacena soluciones aceptadas
ky=0;
kz=0;
k=0;
Error = zeros(Ny-v, Nx-v);
for c1 = 1:Ny-v %menos tamaño ventana para no explorar fuera dominio obs
    for c2 = 1:Nx-v
        Vfx = fx(c1:c1+v-1,c2:c2+v-1); %ventana datos obs
        Vfy = fy(c1:c1+v-1,c2:c2+v-1);
        Vfz = fz(c1:c1+v-1,c2:c2+v-1);
        VX = X(c1:c1+v-1,c2:c2+v-1);
        VY = Y(c1:c1+v-1,c2:c2+v-1);
        Vdxfx = dxfx(c1:c1+v-1,c2:c2+v-1);
        Vdyfx = dyfx(c1:c1+v-1,c2:c2+v-1);
        Vdzfx = dzfx(c1:c1+v-1,c2:c2+v-1);
        Vdxfy = dxfy(c1:c1+v-1,c2:c2+v-1);
        Vdyfy = dyfy(c1:c1+v-1,c2:c2+v-1);
        Vdzfy = dzfy(c1:c1+v-1,c2:c2+v-1);
        Vdxfz = dxfz(c1:c1+v-1,c2:c2+v-1);
        Vdyfz = dyfz(c1:c1+v-1,c2:c2+v-1);
        Vdzfz = dzfz(c1:c1+v-1,c2:c2+v-1);
        %Cambiar de matriz a vector columna
        Vfx = Vfx(:);
        Vfy = Vfy(:);
        Vfz = Vfz(:);
        VX = VX(:);
        VY = VY(:);
        Vdxfx = Vdxfx(:);
        Vdyfx = Vdyfx(:);
        Vdzfx = Vdzfx(:);
        Vdxfy = Vdxfy(:);
        Vdyfy = Vdyfy(:);
        Vdzfy = Vdzfy(:);
        Vdxfz = Vdxfz(:);
        Vdyfz = Vdyfz(:);
        Vdzfz = Vdzfz(:);
        IE = IE(:);
        %FORMANDO SISTEMA DE EC LINEALES SOBREDETERMINADO
        %G:
        GX = [Vdxfx Vdyfx Vdzfx IE];
        GY = [Vdxfy Vdyfy Vdzfy IE];
        GZ = [Vdxfz Vdyfz Vdzfz IE];
        
        GXT = [Vdxfx Vdyfx Vdzfx IE IE*0 IE*0];
        GYT = [Vdxfy Vdyfy Vdzfy IE*0 IE IE*0];
        GZT = [Vdxfz Vdyfz Vdzfz IE*0 IE*0 IE];
        G = [GXT;GYT;GZT];
        %d:
        dfx = VX.*Vdxfx+VY.*Vdyfx+Vfx.*IE;
        dfy = VX.*Vdxfy+VY.*Vdyfy+Vfy.*IE;
        dfz = VX.*Vdxfz+VY.*Vdyfz+Vfz.*IE;
        d = [dfx;dfy;dfz];
        %RESOLVIENDO POR MINIMOS CUADRADOS AMORTIGUADOS d=G*m 
        %m=(G'G+I*e)^(-1)*(G'*d)
        %A=(G'G+I*e)^(-1) e=número pequeño para que no se indetermine
        %inversa
        Ax = inv(GX'*GX+eye(4)*10^(-10)); %comando eye=Matriz Identidad 
        Ay = inv(GY'*GY+eye(4)*10^(-10));
        Az = inv(GZ'*GZ+eye(4)*10^(-10));
        A = inv(G'*G+eye(6)*10^(-10));
        mx = Ax*(GX'*dfx);
        my = Ay*(GY'*dfy);
        mz = Az*(GZ'*dfz);
        m = A*(G'*d);
        
        %APLICANDO CRITERIO THOMPSON (ESTADÍSTICA) Error= d-G*m
        ex = dfx - GX*mx; %vector de errores
        ey = dfy - GY*my;
        ez = dfz - GZ*mz;
        e = d - G*m;
        %error cuadratico E=sum(e.^2)
        Ex = ex'*ex; %suma de los elementos del vector errores elevado cuad
        Ey = ey'*ey;
        Ez = ez'*ez;
        E = e'*e;
        Error(c1, c2) = E;
        %A(3,3): desviación estándar de incognita z0
        sigmax = sqrt((Ax(3,3))*Ex/((v^2)-4)); %4:numero incognitas
        sigmay = sqrt((Ay(3,3))*Ey/((v^2)-4));
        sigmaz = sqrt((Az(3,3))*Ez/((v^2)-4));
        sigma = sqrt((A(3,3))*E/((v^2)-6));
        numx = mx(3)/abs(nu*sigmax); %m(3):posición en z, solucion z0
        numy = my(3)/abs(nu*sigmay);
        numz = mz(3)/abs(nu*sigmaz);
        num = m(3)/abs(nu*sigma);
        if numx>Tol
            %se acepta solución
            kx=kx+1;
            Mx(kx,:)=mx;
        end
        if numy>Tol
            %se acepta solución
            ky=ky+1;
            My(ky,:)=my;
        end
        if numz>Tol
            %se acepta solución
            ky=ky+1;
            Mz(ky,:)=mz;
        end
        if num>Tol
            %se acepta solución
            k=k+1;
            M(k,:)=m;
        end
        figure(2)
        mesh(X,Y,fz)
        hold on
        %surf(VX,VY,Mz)
        view(2)
        hold off
        %pause
    end
end
%M matriz que almacena soluciones que cumplen criterio Thompson
%GRAFICAS:
sxfx = Mx(:,1)*dx;
syfx = Mx(:,2)*dy;
szfx = Mx(:,3)*(dx+dy)*0.5*(-1); %-1 por eje de referencia(coordenadas fuente)
sxfy = My(:,1)*dx;
syfy = My(:,2)*dy;
szfy = My(:,3)*(dx+dy)*0.5*(-1);
sxfz = Mz(:,1)*dx;
syfz = Mz(:,2)*dy;
szfz = Mz(:,3)*(dx+dy)*0.5*(-1);
sx = M(:,1)*dx;
sy = M(:,2)*dy;
sz = M(:,3)*(dx+dy)*0.5*(-1);
%DIMENSIONANDO DOMINIO DE OBSERVACIÓN
X = X*dx;
Y = Y*dy;
figure(3)
contour(X,Y,fx);
hold on
plot3(sxfx,syfx,szfx,'*');
zlim([-1 0])
%title(['Deconvolución de Euler convencional en GX con Ventanas= ',num2str(v),' y Tolerancia=',num2str(Tol)]);
xlabel('Distancia x[km]');
ylabel('Distancia y[km]');
zlabel('Profundidad');
%figure(3)
hold on
contour(X,Y,fy);
hold on
plot3(sxfy,syfy,szfy,'*');
zlim([-1 0])
%title(['Deconvolución de Euler convencional en GY con Ventanas= ',num2str(v),' y Tolerancia=',num2str(Tol)]);
xlabel('Distancia x[km]');
ylabel('Distancia y[km]');
zlabel('Profundidad');
%figure(4)
hold on
contour(X,Y,fz);
hold on
plot3(sxfz,syfz,szfz,'*');
zlim([-1 0])
%title(['Deconvolución de Euler convencional en GZ con Ventanas= ',num2str(v),' y Tolerancia=',num2str(Tol)]);
xlabel('Distancia x[km]');
ylabel('Distancia y[km]');
zlabel('Profundidad'); 
title(['Deconvolución de Euler Convencional con Ventanas= ',num2str(v),' y Tolerancia=',num2str(Tol)]);
legend('gx','GX','gy','GY','gz','GZ');

figure(5)
hold on
contour(X,Y,fz);
hold on
contour(X,Y,fx)
hold on
contour(X,Y,fy)
hold on
plot3(sx,sy,sz,'*');
zlim([-1 0])
%title('Deconvolución de Euler extendida');
title(['Deconvolución de Euler extendida con Ventanas =',num2str(v),' y Tolerancia=',num2str(Tol)]);
xlabel('Distancia x[km]');
ylabel('Distancia y[km]');
zlabel('Profundidad');
%aproxima a fuentes que causan anomalías, describe estructuras
%ejes km

%tol, indice estructural y tamaño ventana
% Graficar el error
figure;
surf(Error);
title('Error en la deconvolución de Euler');
xlabel('Índice en Y');
ylabel('Índice en X');
zlabel('Error');