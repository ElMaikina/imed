%
% Funcion que genera la reconstruccion con filtro, sin filtro
% y con la transformada de Fourier al mismo tiempo, guardando
% el valor de cada reconstruccion en una variable diferente.
%
function [recon_sin_filtro, recon_filtrada, recon_fourier] = retroproyeccion_MS_SS(sinograma, angulos)
    % Convierte los valores de la matriz en floats dobles
    sinograma = double(sinograma);
    % Obtiene las dimensiones del sinograma de entrada
    [num_detectores, ~] = size(sinograma);
    num_angulos = size(sinograma, 2);
    angulos = linspace(0, 180, num_angulos);
    % Tamaño de imagen reconstruida (debe ser cuadrada)
    N = num_detectores;
    % Inicializar imágenes
    recon_sin_filtro = zeros(N, N);
    recon_filtrada   = zeros(N, N);
    % Centro de la imagen
    centro = floor(N/2);
    
    % 1. Back-Projection (sin filtro)
    for i = 1:num_angulos
        theta = angulos(i) * pi / 180;
        proyeccion = sinograma(:, i);
        for x = 1:N
            for y = 1:N
                % Coordenadas centradas
                x_c = x - centro;
                y_c = y - centro;
                % Coordenada proyectada
                t = round(x_c * cos(theta) + y_c * sin(theta)) + centro;
                % Validar rango
                if t >= 1 && t <= num_detectores
                    recon_sin_filtro(x,y) = recon_sin_filtro(x,y) + proyeccion(t);
                end
            end
        end
    end
    % Normalización
    recon_sin_filtro = recon_sin_filtro / num_angulos;
    
    % 2.1. Filtrado (usando rampa en Fourier)
    sinograma_filtrado = zeros(size(sinograma));
    for i = 1:num_angulos
        proyeccion = sinograma(:, i);
        % FFT de la proyección
        P = fft(proyeccion);
        % Crear filtro rampa
        Nf = length(P);
        filtro = zeros(Nf,1);
        for k = 1:Nf
            frecuencia = abs(k - Nf/2);
            filtro(k) = frecuencia;
        end
        % Aplicar filtro
        P_filtrado = P .* filtro;
        % Transformada inversa
        proyeccion_filtrada = real(ifft(P_filtrado));
        sinograma_filtrado(:, i) = proyeccion_filtrada;
    end
    
    % 2.2. Back-Projection (con filtro)
    for i = 1:num_angulos
        theta = angulos(i) * pi / 180;
        proyeccion = sinograma_filtrado(:, i);
        for x = 1:N
            for y = 1:N
                x_c = x - centro;
                y_c = y - centro;
                t = round(x_c * cos(theta) + y_c * sin(theta)) + centro;
                if t >= 1 && t <= num_detectores
                    recon_filtrada(x,y) = recon_filtrada(x,y) + proyeccion(t);
                end
            end
        end
    end
    recon_filtrada = recon_filtrada / num_angulos;

    % 3. Reconstruccion por Fourier
    espacio_k = zeros(N, N);
    for i = 1:num_angulos
        theta = angulos(i) * pi / 180;
        proyeccion = sinograma_filtrado(:, i);
        % FFT de la proyección
        P = fft(proyeccion);
        for k = 1:N
            % Coordenadas en frecuencia
            u = round(centro + (k - centro) * cos(theta));
            v = round(centro + (k - centro) * sin(theta));
            if u >= 1 && u <= N && v >= 1 && v <= N
                espacio_k(u,v) = P(k);
            end
        end
    end
    % Transformada inversa 2D
    recon_fourier = real(ifft2(espacio_k));
end