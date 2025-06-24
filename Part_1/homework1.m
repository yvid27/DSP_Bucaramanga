clearvars
close all
clc
% Maycol's code

% Parametros
f0 = 1e6; % Frecuencia central del pulso (1 MHz)
w0 = 2*pi*f0; % Frecuencia angular (rad/s)

c0 = 1500; % Velocidad del sonido en agua (m/s)
lambda = c0/f0; % Longitud de onda correspondiente (0.0015 m)
rho = 1000;

% Dimensiones del dominio
Lx = 0.1; % 10 cm
Ly = 0.1; % 10 cm
Lz = 0.0015; % 1.5 mm

max_distance = 1.2*sqrt(Lx^2+Ly^2+Lz^2); % distancia diagonal (máxima posible entre fuente y punto)
max_time = max_distance / c0; % tiempo que tarda en propagarse el pulso más lejos posible
f_sample = 10*f0; % frecuencia de muestreo temporal
dt = 1.0/f_sample; % paso temporal
nt = 2*round(max_time /dt/2); % número de muestras temporales
time = (0:nt-10)*dt; % vector temporal
df = 1.0/(nt*dt); % resolución en frecuencia

omega = 2*pi*(0:nt/2+1)*df; % vector de frecuencias angulares

scr =[Lx/2, Ly/2, Lz/2]; % fuente en el centro del dominio

pulse = gauspuls(time-5.0/f0,f0); % pulso centrado temporalmente
Pulse = fft(pulse); % FFT del pulso

n = 1000;
% malla con resolución espacial
% /3
dx = lambda/3;
dy = lambda/3;
dz = lambda/3;

X = 0:dx:Lx; Y = 0:dy:Ly; Z=0:dz:dz;  % una sola capa en Z
[x,y,z] = meshgrid(X,Y,Z);
epsilon = 0.000001; % valor pequeño para evitar división por cero cuando r=0
Radius = sqrt((x-scr(1)).^2+(y-scr(2)).^2+(z-scr(3)).^2)+epsilon; % Calcula el radio r desde cada 
% punto del espacio hasta la fuente

%% Compute G
G = zeros(length(X), length(Y), length(Z),length(omega));

for k = 1:length(omega)
    G(:,:,:,k) = Pulse(k)*(1/(4*pi)) * exp(-1j*omega(k)*Radius/c0)./ Radius;
end

%%
% reconstruye una señal compleja simétrica
wave = zeros(length(X),length(Y),length(Z),nt);
wave(:,:,:,1:nt/2)    = G(:,:,:,1:nt/2);
wave(:,:,:,nt/2+2:nt) = conj(G(:,:,:,nt/2:-1:2));
wave = real(ifft(wave,[],4));

%%
clc
wave_db = 20*log10(1e-10+abs(wave)/max(abs(wave(:))));
cc = [-100 -60];

%%
 figure(1);
 for  k = 10:nt
     imagesc(squeeze(wave_db(:,:,end/2,k)))
     caxis(cc)
     pause(0.1)
 end
            return