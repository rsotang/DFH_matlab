%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LA RELATIVA NO ESTÁ BIEN IMPLEMENTADA, AHORA HACE EL % DEL MAXIMO DEL
%ARCHIVO, NO DEL MAXIMO DE LA PRESCRIPCIO
%HAY QUE IMPLEMENTAR O LA LECTURA DE UN RTPLAN O UNA INTRODUCCIÓN DE
%USUARIO




info = dicominfo("nmT6Pa0ScaJ1BpWCn65cfZ004\rt-struct.dcm");
%vol = importar_cabeza();
vol = load("bod.mat");
vol = vol.vol;

%%
rtContours =dicomContours(info);
rtContours.ROIs;
%%
figure
plotContour(rtContours)
ax = gca;

vol_ref = imref3d(size(vol),ax.XLim,ax.YLim,ax.ZLim); %creamos un espacio de referencia en el que se cortan los volumenes (no hay una referencia geometrica entre ambos objetos)
rtMask = createMask(rtContours,5,vol_ref);
close;

%maskDisp = volshow(rtMask); %todavía no se si funciona de todo bien la
%mascara, visualmente es rara 
%%
%AHORA HAY QUE SER CAPAZ DE SACAR INFORMACIÓN RELEVANTE DEL HISTOGRAMA, QUE
%ES DONDE ME HACE FALTA LA INFORMACIÓN DEL VOLUMEN

info_CT = dicominfo("nmT6Pa0ScaJ1BpWCn65cfZ004\CT.nmT6Pa0ScaJ1BpWCn65cfZ004.Image 1.dcm");
x_voxel_size = info_CT.PixelSpacing(1); %en mm
y_voxel_size = info_CT.PixelSpacing(2);
z_voxel_size = info_CT.SliceThickness;
pixel_vol = x_voxel_size+y_voxel_size * z_voxel_size;
total_numel_mask = nnz(rtMask);
total_vol_mask = total_numel_mask * x_voxel_size * y_voxel_size *z_voxel_size / 1000; %hay que verificar estas cuentas

%TAMBIEN HAY QUE MIRAR COMO REPONDERAR LOS VOXELES DE VOLUMEN
%%
%aplicar esa mascara a un archivo de dosis
info_dose=dicominfo("nmT6Pa0ScaJ1BpWCn65cfZ004\rt-dose.dcm");
pixel_data = dicomread("nmT6Pa0ScaJ1BpWCn65cfZ004\rt-dose.dcm");
pixel_data = squeeze(pixel_data); %no recuerdo porque venia con una dimensión mas pero si que habia que apachurrarlos
%IMPORTANTE
%La dosis que viene codificada en los archivos de rtdose como pixeldata no
%tiene un valor absoluto de escala de grises-dosis
%hay una etiqueta DICOM con un factor multiplicativo que pondera el valor
%de pixel y lo convierte en dosis
%tambien tener cuidado con las etiquetas que te dicen si la dosis es FISICA
%o EQUIVALENTE y si es por FRACCION o por PLAN (en nuestro caso fisica y
%plan)
pixel_to_dose_factor = info_dose.DoseGridScaling;
dose = double(pixel_data).* pixel_to_dose_factor;
dose = imresize3(dose,size(rtMask));
%mask = imresize3(rtMask,size(dose));
dose_double = double(dose);
mask_double = double(rtMask);
data_masked_double = dose_double .* mask_double; %si, es así de facil, es la gracia de hacer una mascara binaria
data = data_masked_double; %para programar con más comodidad

%buscamos el valor maximo del array y lo tomamos como maximo de la lista
%la lista estará vacía y en  cuanto encuentre el max tendrá el numero de elementos en 
% cGy del max (i.e si encuentra un max de 57.6 Gy tendremos un array de len=5760

max_value = max(max(max(data)));
rel_data = data .*100/max_value;
%triple bucle anidado para recorrer todo el volumen de la mascara (acotar entre los cortes no vacios) Cada valor de pixel añade una cuenta a todos los huecos de array por debajo o igual que el 
cGy_dose = zeros(1,round(100*max_value));%cGy
percentage_cGy_dose = zeros(1,100);% rel
%representamos el array como un histograma (??????????????)
%voy a hacer un triple bucle recorriendo cada pixel del array enmascarado.
%Cada vez que se lea un valor se sumará +1 a todos los valores anteriores

indices = size(rtMask);

for i=1:indices(1)
    for j=1:indices(2)
        for k=1:indices(3)
          dose_pixel_value = data(i,j,k);
          if dose_pixel_value == 0
              continue
          else          
          for q=1:round(100*dose_pixel_value)
          cGy_dose(1,q)= cGy_dose(1,q) + 1;
          
          end
          end
        end
    end
end
%%%%%%%%%%REL%%%%%%%%%%%%%%
for i=1:indices(1)
    for j=1:indices(2)
        for k=1:indices(3)
          dose_pixel_value = rel_data(i,j,k);
          if dose_pixel_value == 0
              continue
          else          
          for q=1:round(dose_pixel_value)
          percentage_cGy_dose(1,q)= percentage_cGy_dose(1,q) + 1;
          
          end
          end
        end
    end
end
%renormalizamos
percentage_cGy_dose_100 = percentage_cGy_dose .* 100/percentage_cGy_dose(1,1);
percentage_cGy_dose_cm3 = percentage_cGy_dose .* pixel_vol /1000;

% Crear el eje x, que representaría la dosis
dosis = 1:numel(cGy_dose);  

categories_dose_cm3= cGy_dose .* pixel_vol / 1000;
categories_dose_rel= cGy_dose .* 100/cGy_dose(1,1);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AHORA LA LIAMOS PARDA CON EL SPECT
%%
%importamos SPECT en un array 3d
%[spect,info_spect]= importar_stack('SPECT');
spect = load("D:\Polyspace\R2020b\=SCRIPTS=\mascara rt_struct\spect.mat");
%spect = load("H:\Polyspace\R2020b\=SCRIPTS=\mascara rt_struct\spect.mat");
spect = spect.spect;
%reescalamos al tamaño de la estructura
spect = imresize3(spect,size(data));
%enmascaramos tambien el spect con la struct para quedarnos con el pulmon que queremos 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%HAY QUE TENER MUCHO CUIDADO AL HACER EL IMRESIZE, HAY QUE VERIFICAR EL
%%%%%%%%%%%%%%%%%%%%%%%%%%GROSOR DE CORTE Y VER SI INTRODUCIMOS ALGUNA DEFORMACIÓN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
spect_masked = spect .* mask_double;
%normalizamos los valores del SPECT (añadir nota los max normalizarlos a 0.8 si seguimos el paper)
spect_max_value = max(max(max(spect_masked)));
spect_normalized = spect_masked.*1/(0.8*spect_max_value);
%multiplicamos cada pixel de la mascara por su valor de pixel correspondiente del SPECT

functional_data = data .* spect_normalized;
%repetimos la calculador del DHV
%%
max_function = round(100*max(max(max(functional_data))));
cGy_function = zeros(1,max_function);
rel_functional_data = functional_data .* 100/max_function;
percentage_cGy_function = zeros(1,100);

for i=1:indices(1)
    for j=1:indices(2)
        for k=1:indices(3)
            dose_pixel_value = functional_data(i,j,k);
          if dose_pixel_value == 0
              continue
          else          
             for q=1:round(100*dose_pixel_value)
                cGy_function(1,q)= cGy_function(1,q) + 1;
             end
          end
        end
    end
end

%%%%%%%%%%REL%%%%%%%%%%%%%%%%%%

for i=1:indices(1)
    for j=1:indices(2)
        for k=1:indices(3)
          dose_pixel_value = rel_functional_data(i,j,k);
          if dose_pixel_value == 0
              continue
          else          
          for q=1:round(dose_pixel_value)
          percentage_cGy_function(1,q)= percentage_cGy_function(1,q) + 1;
          
          end
          end
        end
    end
end
percentage_cGy_function_100 = percentage_cGy_function .* 100/percentage_cGy_function(1,1);
percentage_cGy_function_cm3 = percentage_cGy_function .* pixel_vol /1000;

nu_dose = 1:numel(cGy_function);  % mismo marcao de referencia que el DVH???
categories_function_cm3 = cGy_function .* pixel_vol /1000;
categories_function_rel = cGy_function .* 100/cGy_function(1,1);
percentage = 1:100;

%% Graficar el DVH


% figure(1);
% %axes_cm3=axes;
% plot(dosis, categories_dose_cm3 , 'LineWidth', 2);%categories dose por que se clasifica en volumen ( en este caso voxeles)
% xlabel('Dosis (cGy)');
% ylabel('Volumen (cm3)');
% title('Histograma Dosis-Volumen Acumulativo (DVH)');
% grid on;
% 
% figure(2);
% %axes_voxel= axes;
% plot(percentage, percentage_cGy_dose_cm3 , 'LineWidth', 2);
% xlabel('Dosis (%)');
% ylabel('Volumen (cm3)');
% title('Histograma Dosis-Volumen Acumulativo (DVH)');
% grid on;
% 
% 
% figure(3);
% %axes_rel=axes;
% plot(dosis, categories_dose_rel , 'LineWidth', 2);
% xlabel('Dosis (cGy)');
% ylabel('Volumen (%)');
% title('Histograma Dosis-Volumen Acumulativo (DVH)');
% grid on;
% 
% 
% figure(4);
% %axes_percentage = axes;
% plot(percentage, percentage_cGy_dose_100 , 'LineWidth', 2);
% xlabel('Dosis (%)');
% ylabel('Volumen (%)');
% title('Histograma Dosis-Volumen Acumulativo (DVH)');
% 
% grid on;
% 
% %% Graficar el DFH
% %figure;
% figure(1);
% hold on;
% plot(nu_dose, categories_function_cm3, 'LineWidth', 2);
% hold off;
% 
% figure(2);
% hold on;
% plot(percentage, percentage_cGy_function_cm3, 'LineWidth', 2);
% hold off;
% 
% figure(3);
% hold on;
% plot(nu_dose, categories_function_rel, 'LineWidth', 2);
% hold off;
% 
% figure(4);
% hold on
% plot(percentage, percentage_cGy_function_100, 'LineWidth', 2);
% hold off;
%%

CalculadoraValoresHistograma(categories_dose_cm3,categories_dose_rel, percentage_cGy_dose_100,percentage_cGy_dose_cm3,dosis,categories_function_cm3,categories_function_rel,percentage_cGy_function_cm3,percentage_cGy_function_100,nu_dose,percentage)


%ahora necesitamos una funcion que nos saque parametros dosimetricos
%podemos hacer el codigo aqui a lo guarro para que nos de el V20 y luego
%implementarlo como funcion. 
%tambien interesa corregir el numero de pixeles por volumen o sacar en
%relativo. hay que cambiar el valor de los elementos de categories_dose
%tambien interesa hacer esto con un paciente con CT4d porque asi
%implementamos los dos features simultaneamente y además hacemos una
%comprobación clinica real.



%los vectores a usar para extraer valores son categories_cm3,
%categories_rel, percentage_100, percentage_cm3

function CalculadoraValoresHistograma(categories_dose_cm3,categories_dose_rel, percentage_cGy_dose_100,percentage_cGy_dose_cm3,dosis,categories_function_cm3,categories_function_rel,percentage_cGy_function_cm3,percentage_cGy_function_100,nu_dose,percentage)%,percentage)
    % Crear la figura
    hFig = figure('Position', [2100, 300, 1280, 720], 'Name', 'DVH calculator', 'NumberTitle', 'off', 'MenuBar', 'none');
    ax = axes('Parent',hFig,'position',[0.07 0.1 0.7 0.8]);
    p1=plot(ax,dosis,categories_dose_cm3,'LineWidth',2);
    hold on
    p2=plot(ax,nu_dose,categories_function_cm3,'LineWidth',2);
    hold off
    grid on
    xlabel('Dosis (cGy)');
    ylabel('Volumen (cm3)');
    title('Histograma Dosis-Volumen Acumulativo (DVH) + Dosis funcional-Volumen (DFH)')
    
    %%%%%%%%%%%%%%%%%%%ELEMENTOS GRÁFICOS DVH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mensaje de texto
    uicontrol('Style', 'text', 'String', 'Introduce restricción de dosis a calcular:', 'Position', [1050, 650, 200, 20]);
    % Recuadros para introducir datos
    hEdit1 = uicontrol('Style', 'popupmenu', 'Position', [1060, 620, 50, 20], 'String', {'V','D'}, 'Callback', @updateOptions);
    
    hEdit2 = uicontrol('Style', 'edit', 'Position', [1130, 620, 40, 20]);
    
    hEdit3 = uicontrol('Style', 'popupmenu', 'Position', [1190, 620, 50, 20], 'Callback', @updateGraphics);
    hEdit3.String={'Gy','%'};
    
    hEdit4 = uicontrol('Style', 'popupmenu', 'Position', [1190, 500, 50, 20], 'Callback', @updateGraphics);
    hEdit4.String={'cm3','%'};
    % Botón de calcular
    hButton = uicontrol('Style', 'pushbutton', 'String', 'Calcular', 'Position', [1050, 500, 100, 20], 'Callback', @calculate_callback);
    
    %respuesta
    result_button=uicontrol('Style', 'text', 'String', '', 'Position', [1120, 580, 60, 20],'BackgroundColor',[1,1,1]);
    uicontrol('Style', 'text', 'String', 'DVH', 'Position', [1050, 575, 60, 20]);
    fresult_button=uicontrol('Style', 'text', 'String', '', 'Position', [1120, 540, 60, 20],'BackgroundColor',[1,1,1]);
    uicontrol('Style', 'text', 'String', 'DFH', 'Position', [1050, 535, 60, 20]);
    
    
   
    
    %%%%%%%%%%%%%%%%%%%%%%Funciones de interfaz gráfica%%%%%%%%%%%%%%%%%%%%
    
    
    function updateOptions(~,~)
        switch hEdit1.String{hEdit1.Value}
            case 'V'
                hEdit3.String = {'Gy','%'};
                hEdit4.String = {'cm3','%'};
            case 'D'
                hEdit3.String = {'cm3','%'};
                hEdit4.String = {'Gy','%'};
        end
    end

    function updateGraphics(~,~)
      switch hEdit1.String{hEdit1.Value}
            case 'V' %puede tener 
                if strcmp(hEdit4.String{hEdit4.Value},'cm3') && strcmp(hEdit3.String{hEdit3.Value},'Gy')
                    p1=plot(ax,dosis,categories_dose_cm3,'LineWidth',2);
                    hold on
                    p2=plot(ax,nu_dose,categories_function_cm3,'LineWidth',2);
                    hold off
                    grid on
                    xlabel('Dosis (cGy)');
                    ylabel('Volumen (cm3)');
                    title('Histograma Dosis-Volumen Acumulativo (DVH) + Dosis funcional-Volumen (DFH)')
                elseif strcmp(hEdit4.String{hEdit4.Value},'%') && strcmp(hEdit3.String{hEdit3.Value},'Gy')
                    p1=plot(ax,dosis,categories_dose_rel,'LineWidth',2);
                    hold on
                    p2=plot(ax,nu_dose,categories_function_rel,'LineWidth',2);
                    hold off
                    grid on
                    xlabel('Dosis (cGy)');
                    ylabel('Volumen (%)');
                    title('Histograma Dosis-Volumen Acumulativo (DVH) + Dosis funcional-Volumen (DFH)')
                elseif strcmp(hEdit4.String{hEdit4.Value},'cm3') && strcmp(hEdit3.String{hEdit3.Value},'%')
                    p1=plot(ax,percentage, percentage_cGy_dose_cm3,'LineWidth',2);
                    hold on
                    p2=plot(ax,percentage,percentage_cGy_function_100,'LineWidth',2);
                    hold off
                    grid on
                    xlabel('Dosis (%)');
                    ylabel('Volumen (cm3)');
                    title('Histograma Dosis-Volumen Acumulativo (DVH) + Dosis funcional-Volumen (DFH)')
                elseif strcmp(hEdit4.String{hEdit4.Value},'%') && strcmp(hEdit3.String{hEdit3.Value},'%')
                    p1=plot(ax,percentage,percentage_cGy_dose_100,'LineWidth',2);
                    hold on
                    p2=plot(ax,percentage,percentage_cGy_function_100,'LineWidth',2);
                    hold off
                    grid on
                    xlabel('Dosis (%)');
                    ylabel('Volumen (%)');
                    title('Histograma Dosis-Volumen Acumulativo (DVH) + Dosis funcional-Volumen (DFH)')
                end
            case 'D' % puede %%, %Gy, cm3Gy, cm3%
                 if strcmp(hEdit3.String{hEdit3.Value},'cm3') && strcmp(hEdit4.String{hEdit4.Value},'Gy')
                    p1=plot(ax,dosis,categories_dose_cm3,'LineWidth',2);
                    hold on
                    p2=plot(ax,nu_dose,categories_function_cm3,'LineWidth',2);
                    hold off
                    grid on
                    xlabel('Dosis (cGy)');
                    ylabel('Volumen (cm3)');
                    title('Histograma Dosis-Volumen Acumulativo (DVH) + Dosis funcional-Volumen (DFH)')
                elseif strcmp(hEdit3.String{hEdit3.Value},'%') && strcmp(hEdit4.String{hEdit4.Value},'Gy')
                    p1=plot(ax,dosis,categories_dose_rel,'LineWidth',2);
                    hold on
                    p2=plot(ax,nu_dose,categories_function_rel,'LineWidth',2);
                    hold off
                    grid on
                    xlabel('Dosis (cGy)');
                    ylabel('Volumen (%)');
                    title('Histograma Dosis-Volumen Acumulativo (DVH) + Dosis funcional-Volumen (DFH)')
                elseif strcmp(hEdit3.String{hEdit3.Value},'cm3') && strcmp(hEdit4.String{hEdit4.Value},'%')
                    p1=plot(ax,percentage, percentage_cGy_dose_cm3,'LineWidth',2);
                    hold on
                    p2=plot(ax,percentage,percentage_cGy_function_100,'LineWidth',2);
                    hold off
                    grid on
                    xlabel('Dosis (%)');
                    ylabel('Volumen (cm3)');
                    title('Histograma Dosis-Volumen Acumulativo (DVH) + Dosis funcional-Volumen (DFH)')
                elseif strcmp(hEdit3.String{hEdit3.Value},'%') && strcmp(hEdit4.String{hEdit4.Value},'%')
                    p1=plot(ax,percentage,percentage_cGy_dose_100,'LineWidth',2);
                    hold on
                    p2=plot(ax,percentage,percentage_cGy_function_100,'LineWidth',2);
                    hold off
                    grid on
                    xlabel('Dosis (%)');
                    ylabel('Volumen (%)');
                    title('Histograma Dosis-Volumen Acumulativo (DVH) + Dosis funcional-Volumen (DFH)')
                end
        end
    end



    %%%%%%%%%%%%%%%% Función de callback para el botón DVH%%%%%%%%%%%%%%%%%
    
    function calculate_callback(~,~)
    % Recuperar los valores introducidos por el usuario
    variable = hEdit1.String{hEdit1.Value};
    val = str2num(hEdit2.String);
    units = hEdit3.String{hEdit3.Value};
    units_result = hEdit4.String{hEdit4.Value};
    
    result = 0;  % inicializar resultado
    fresult = 0; % inicializar fresultado
    try
    %si V y Gy%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(variable, 'V') && strcmp(units, 'Gy')
      val = 100 * val; %pasamos a cGy     
        if strcmp(units_result,'cm3')           
            result = round(categories_dose_cm3(val),2);          
           
            fresult = round(categories_function_cm3(val),2);
            
        elseif strcmp(units_result,'%')
            result = round(categories_dose_rel(val),2);
            fresult = round(categories_function_rel(val),2);
        end
      
    %si V y %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif strcmp(variable, 'V') && strcmp(units, '%')
        
        if strcmp(units_result,'cm3')
            result = round(percentage_cGy_dose_cm3(val),2);
            fresult = round(percentage_cGy_function_cm3(val),2);
        elseif strcmp(units_result,'%')
             result = round(percentage_cGy_dose_100(val),2);
             fresult = round(percentage_cGy_function_100(val),2);
        end
        
    %si D y %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif strcmp(variable, 'D') && strcmp(units, '%')
        if strcmp(units_result,'Gy')
           result = round(HistogramaInterpolado(percentage_cGy_dose_cm3,val)/100,2);
           fresult = round(HistogramaInterpolado(percentage_cGy_function_cm3,val)/100,2);
        elseif strcmp(units_result,'%')
           result = round(HistogramaInterpolado(percentage_cGy_dose_100,val),2);
           fresult = round(HistogramaInterpolado(percentage_cGy_function_100,val),2);
        end

    %si D y cm3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif strcmp(variable, 'D') && strcmp(units, 'cm3')
        if strcmp(units_result,'Gy')
            result = round(HistogramaInterpolado(categories_dose_cm3,val)/100,2);
            fresult = round(HistogramaInterpolado(categories_function_cm3,val)/100,2);
        elseif strcmp(units_result,'%')
            result = round(HistogramaInterpolado(categories_dose_rel,val),2);
            fresult = round(HistogramaInterpolado(categories_function_rel,val),2);
        end

    end
    
    % Esta parte supone que defines 'result' en uno de los bloques anteriores.
    result_button.String = num2str(result);
    fresult_button.String = num2str(fresult);
    
    catch ME
        if strcmp(ME.identifier, 'MATLAB:badsubscript') || strcmp(ME.identifier, 'MATLAB:badsubstruct')
            % Si es un error de índice, establece el resultado a 0.
            result_button.String = '0';
            fresult_button.String = '0';
        else
            % Para cualquier otro error, muestra el mensaje de error.
            disp(['Error: ' ME.message]);
        end
    end
    
end

function result = HistogramaInterpolado(array, valorBuscado)
    
    % Encontrar la primera posición en el array cuyo valor sea menor que el valor buscado
    idx_inferior = find(array <= valorBuscado, 1, 'first');
    
    % Si el valor buscado es mayor que el valor máximo en el array, 
    % o es menor que el valor mínimo, no podemos hacer la interpolación.
    if isempty(idx_inferior) || idx_inferior == 1 || idx_inferior == length(array)
        error('Valor buscado fuera de rango para interpolación');
    end
    
    % El índice del valor superior es simplemente el índice inferior menos 1
    idx_superior = idx_inferior - 1;
    
    % Valores para la interpolación
    y1 = array(idx_superior);
    y2 = array(idx_inferior);
    x1 = idx_superior;
    x2 = idx_inferior;
    
    % Interpolación lineal
    result = x1 + (x2 - x1) * (valorBuscado - y1) / (y2 - y1);
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------TO--DO------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-DESCARGAR PACIENTE CON CT4D Y CAMBIARLO POR EL ACTUAL
%
%-IMPLEMENTAR TAMAÑO DE VECTORES EN FUNCION DE DOSIS MAXIMA PRESCRITA EN
%PLAN (ETIQUETA DICOM) PARA RELATIVA
