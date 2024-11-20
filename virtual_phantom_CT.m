%Voy a generar una serie de maniquís de imagenes sinteticas  para la
%validación del software. Basicamente vamos a generar imagenes de CT, SPECT
%y CT4D sintéticas que cumplan las caracteristicas que necistamos para cada
%imagen y permitan ver las fallas de la implementación del codigo. Para
%ello vamos a generar esferas de distinta resolución y con distintos
%diametros. Es decir, el CT será un bloque 512x512x25 en el que el
%extertior de la esfera tendrá valor -1000 el contornto exterior de la
%esfera puede tener valor 0 y el interior -500. Para el spect dibujamos
%solo media esfera, y para el CT4D una esfera mas grande con el volumen más
%fino y viceversa para los dos bloques del volumen respiratorio. La idea
%también es poder volcar esos maniquis en objetos DICOM para importar a
%Eclipse y hacer la validación cruzada.

%haciendo una esfera hueca con corteza para simular el pulmón

%%
% Parámetros del volumen y de la esfera
volSize = [512, 512, 512];   % Tamaño del volumen
outerRadius = 125;  % Radio exterior de la corteza esférica
innerRadius = outerRadius * 0.8; % Radio interior (define el grosor de la corteza)

% Generar el volumen
volumeCT = generateSphereCT(volSize, outerRadius, innerRadius);

% Mostrar el volumen 3D
figure;
volshow(volumeCT);  % Mostrar el volumen en 3D (requerido Matlab R2020b o posterior)

% Mostrar el corte central
showCentralSlice(volumeCT);  % Llamar a la función para visualizar el corte central
%%
% Función principal que genera el volumen con una corteza esférica y valores específicos
function volume = generateSphereCT(volSize, outerRadius, innerRadius)
    % Inicializar el volumen lleno de -1000 (aire)
    volume = -1000 * ones(volSize);
    center = volSize / 2;  % Centro de la esfera

    % Rellenar el volumen con la corteza esférica y el interior
    for z = 1:volSize(3)
        for y = 1:volSize(2)
            for x = 1:volSize(1)
                % Distancia euclidiana al centro de la esfera
                dist = sqrt((x - center(1))^2 + (y - center(2))^2 + (z - center(3))^2);

                % Si está dentro del radio exterior pero fuera del radio interior
                if dist <= outerRadius && dist >= innerRadius
                    volume(x, y, z) = 0;  % Corteza esférica con valor 0
                elseif dist < innerRadius
                    volume(x, y, z) = -500;  % Interior de la esfera con valor -500
                end
            end
        end
    end
end

% Función que muestra un corte central del volumen en el plano axial
function showCentralSlice(volume)
    centerSlice = round(size(volume, 3) / 2);  % Índice del corte central en el eje Z
    figure;
    imagesc(volume(:, :, centerSlice));  % Mostrar la imagen del corte central
    colormap('gray');  % Colormap de grises, típico en imágenes CT
    colorbar;  % Barra de colores
    title('Corte central de la esfera (Plano Axial)');
    axis equal;
end



