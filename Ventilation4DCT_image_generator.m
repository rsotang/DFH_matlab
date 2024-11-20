%4DCT_ventilation_image_generator
%vamos a replicar el paper de Dynamic ventilation imaging from four-dimensional
%computed tomography doi:10.1088/0031-9155/51/4/002 en MATLAB
%%
%El registro
%Pruebas de registro rígido
v1=load('Validation_CT0'); % Esfera 1
v2=load('Validation_CT100'); % Esfera 2

volumeCT0= v1.volumeCT0;
volumeCT100= v2.volumeCT100;

%vamos a hacer downsampling de los volumenes para facilitar los cálculos
%con matrices más pequeñas.

% Crear objetos de volúmenes
fixedVolume = imref3d(size(volumeCT0));
movingVolume = imref3d(size(volumeCT100));

% Registro inicial usando el algoritmo predeterminado (afín o deformable)
[optimizer, metric] = imregconfig('multimodal');
movingReg = imregister(volumeCT100, movingVolume, volumeCT0, fixedVolume, 'affine', optimizer, metric);

% Mostrar los volúmenes
volshow(movingReg);
%con estos parametros el registro rígido falla, se queda sin memoria por el
%número de iteraciones.
%%%%%%%%%%%
%%
%vamos a probar el registro defomrable con Optical Flow antes de ponernos a
%ajustar parametros
iterations = 10;  % Número de iteraciones para la convergencia
alpha = 1.0;      % Peso regularizador
[dimX, dimY, dimZ] = size(volumeCT0);

% Inicializar desplazamientos en las 3 dimensiones
Ux = zeros(dimX, dimY, dimZ);
Uy = zeros(dimX, dimY, dimZ);
Uz = zeros(dimX, dimY, dimZ);

% Optical Flow basado en diferencias locales
for iter = 1:iterations
    % Calcular gradientes espaciales de volume2
    [Ix, Iy, Iz] = gradient(volumeCT100);
    
    % Diferencia de intensidad entre imágenes
    It = volumeCT100 - volumeCT0;
    
    % Regularización del desplazamiento
    Ux_smooth = imgaussfilt3(Ux, alpha);
    Uy_smooth = imgaussfilt3(Uy, alpha);
    Uz_smooth = imgaussfilt3(Uz, alpha);
    
    % Actualización de desplazamientos (ecuaciones del flujo óptico)
    num = -(Ix .* Ux_smooth + Iy .* Uy_smooth + Iz .* Uz_smooth + It);
    denom = Ix.^2 + Iy.^2 + Iz.^2 + alpha;
    Ux = Ux + num .* Ix ./ denom;
    Uy = Uy + num .* Iy ./ denom;
    Uz = Uz + num .* Iz ./ denom;
end

% Aplicar el campo de deformación a volume2
[x, y, z] = ndgrid(1:dimX, 1:dimY, 1:dimZ);
deformedVolume = interp3(volumeCT100, x + Ux, y + Uy, z + Uz, 'linear', 0);

% Visualizar resultado
volshow(deformedVolume);

%%
%cálculo voxel a voxel y promediado a 9x9
