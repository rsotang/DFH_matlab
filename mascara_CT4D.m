%&Vamos a hacer lo mismo que para el spect pero ahora para una imagen de
%ventilacion de tac, primero vamos a adaptar la parte de importacion de
%imagenes para los dos CT
info0 = dicominfo("REGISTRO\rt-struct.dcm");
info_CT = dicominfo("nmT6Pa0ScaJ1BpWCn65cfZ004\CT.nmT6Pa0ScaJ1BpWCn65cfZ004.Image 1.dcm");
x_voxel_size = info_CT.PixelSpacing(1); %en mm
y_voxel_size = info_CT.PixelSpacing(2);
z_voxel_size = info_CT.SliceThickness;
pixel_vol = x_voxel_size+y_voxel_size * z_voxel_size;
%info50 = dicominfo("Lungs_50\rt-struct.dcm");
%vol = importar_stack('stack');
vol0 = load("CT_0.mat");
vol0 = vol0.vol0;
vol50 = load("CT_50.mat");
vol50 = vol50.vol50;

%%
%importamos las estructuras de los dos CT
rtContours0 =dicomContours(info0);
rtContours0.ROIs;
%rtContours50 =dicomContours(info50);
%rtContours50.ROIs;
%%
%Creamos la mascara para el primer CT del lung_L
figure
plotContour(rtContours0)
ax = gca;

vol_ref = imref3d(size(vol0),ax.XLim,ax.YLim,ax.ZLim); %creamos un espacio de referencia en el que se cortan los volumenes (no hay una referencia geometrica entre ambos objetos)
rtMask_lungL = createMask(rtContours0,2,vol_ref);
close;
%%
%Creamos la mascara para el lung_R
figure
plotContour(rtContours0)
ax = gca;

vol_ref = imref3d(size(vol50),ax.XLim,ax.YLim,ax.ZLim); %creamos un espacio de referencia en el que se cortan los volumenes (no hay una referencia geometrica entre ambos objetos)
rtMask_lungR = createMask(rtContours0,3,vol_ref);
close;

%%
%Recortamos lung_L y lung_R en CT_0 y CT_50
lung_L_0 = vol0 .*rtMask_lungL;
lung_L_50= vol50 .* rtMask_lungL;
lung_R_0= vol0 .*rtMask_lungR;
lung_R_50= vol50 .* rtMask_lungR;
%%
%Generamos imagen de ventilación segun el paper de Houston

vent_lung_L=zeros(size(lung_L_0));
indices = size(vent_lung_L);
for i=1:indices(1)
    for j=1:indices(2)
        for k=1:indices(3)           
            vent_lung_L(i,j,k) = 1000*(lung_L_0(i,j,k)-lung_L_50(i,j,k))/((lung_L_50(i,j,k)*(1000+lung_L_0(i,j,k))));
                       
        end
    end
end
vent_lung_L(isnan(vent_lung_L)) = 0;
vent_lung_L(vent_lung_L==Inf)=1000;
vent_lung_L(vent_lung_L==-Inf)=-1000;
nu_vent_lung_L=replaceExtremes(vent_lung_L,999);
%si hacemos 1/0 nos da Inf y si hacemos 0/0 nos da NaN
%%NO TENGO NI IDEA DE SI AFECTA O DE DONDE SALEN EL INF Y EL NAN
vent_lung_R=zeros(size(lung_R_0));

for i=1:indices(1)
    for j=1:indices(2)
        for k=1:indices(3)           
            vent_lung_R(i,j,k) = 1000*(lung_R_0(i,j,k)-lung_R_50(i,j,k))/(lung_R_50(i,j,k)*(1000+lung_R_0(i,j,k)));
              
        end
    end
end

vent_lung_R(isnan(vent_lung_R)) = 0;
vent_lung_R(vent_lung_L==Inf)=1000;
vent_lung_R(vent_lung_L==-Inf)=-1000;
nu_vent_lung_R=replaceExtremes(vent_lung_R,999);

global vent
vent = vent_lung_L + vent_lung_R;

%%
hFig = figure('Name', 'Navegador de Array 3D');


% Mostrar el primer corte
current_slice = 1;
h=imshow(vent(:, :, current_slice));
colormap(jet)
colorbar
% Crear el slider
hSlider = uicontrol('Parent',hFig,'Style', 'slider', 'Min', 1, 'Max', size(vent, 3), 'Value', current_slice, ...
    'Units', 'normalized', 'Position', [0.1 0.01 0.8 0.05]);
hSlider.Callback = @(source,ed)set(h,'CData',vent(:,:,round(source.Value)));
%ni puta idea de como funciona esto (copiado de stackoverflow)


%%
%%ESTO PARA CUANDO TENGAMOS HECHA LA IMAGEN DE VENTILACIÓN%%%%%%%%%%%%%
%Importamos el archivo de dosis para más adelante
info_dose=dicominfo("Ventilation\rt-dose.dcm");
pixel_data = dicomread("Ventilation\rt-dose.dcm");
pixel_data = squeeze(pixel_data);
pixel_to_dose_factor = info_dose.DoseGridScaling;
dose = double(pixel_data).* pixel_to_dose_factor;
dose = imresize3(dose,size(rtMask_lungL));



function B = replaceExtremes(A, threshold)
    % A: matriz 3D
    % threshold: umbral para identificar valores extremos (e.g., 1000)
    
    B = A; % matriz de salida
    [rows, cols, pages] = size(A);
    
    for x = 2:rows-1
        for y = 2:cols-1
            for z = 2:pages-1
                % Si el valor es extremo
                if abs(A(x, y, z)) > threshold
                    % Extraer el cubo 3x3x3 alrededor del punto extremo
                    cube = A(x-1:x+1, y-1:y+1, z-1:z+1);
                    % Reemplazar el valor extremo con el promedio del cubo
                    B(x, y, z) = mean(cube(:));
                end
            end
        end
    end
end
%%%%%%%TO DO%%%%%%
%IMPLEMENTAR EN EL CALCULO DOSIMÉTRICO
%PENSAR EN UN MAPA DE COLORES NUEVO
%ENTENDER COMO FUNCIONAN LOS DICOM DE REGISTROS
%INTENTAR HACER EL REGISTRO FUERA, ES DECIR, HACER EL REGISTRO EN ECLIPSE Y
%DESCARGAR LOS CTS SIN REGISTRAR Y APLICAR EL REGISTRO EN MATLAB
%INTENTAR GENERAR REGISTROS DE FORMA AUTONOMA