%
% Reconstruye una imagen con ruido usando la Transformada Discreta de
% Fourier (TDF). Primero calcula el TDF, luego crea una mascara de filtrado
% la cual usa para quitar el ruido y finalmente regenera la imagen usando
% la TDF inversa.
%
function imagen_reconstruida = reconstruccion(imagen)
    % Convertir a double para precisión
    imagen = double(imagen);
    % Obtener dimensiones
    [Nx, Ny] = size(imagen);
    % Calcular TDF
    espectro = zeros(Nx, Ny);
    for u = 1:Nx
        for v = 1:Ny
            suma = 0;
            for x = 1:Nx
                for y = 1:Ny
                    exponente = -2*pi*1i * ( ...
                        (u-1)*(x-1)/Nx + (v-1)*(y-1)/Ny );
                    suma = suma + imagen(x,y) * exp(exponente);
                end
            end
            espectro(u,v) = suma;
        end
    end
    % Crear máscara de filtrado
    mascara = ones(Nx, Ny);
    centro_x = round(Nx/2);
    centro_y = round(Ny/2);
    % Eliminamos puntos equiespaciados en el espectro
    separacion = 10;
    for u = 1:separacion:Nx
        for v = 1:separacion:Ny
            mascara(u,v) = 0;
        end
    end
    % Eliminar líneas diagonales
    for u = 1:Nx
        for v = 1:Ny
            % condición de diagonal
            if abs((u - centro_x) - (v - centro_y)) < 3
                mascara(u,v) = 0;
            end
        end
    end
    % Se aplica el filtro
    espectro_filtrado = espectro .* mascara;
    % Calcular TDF inversa (ahora sin el ruido original)
    imagen_reconstruida = zeros(Nx, Ny);
    for x = 1:Nx
        for y = 1:Ny
            suma = 0;
            for u = 1:Nx
                for v = 1:Ny
                    exponente = 2*pi*1i * ( ...
                        (u-1)*(x-1)/Nx + (v-1)*(y-1)/Ny );
                    suma = suma + espectro_filtrado(u,v) * exp(exponente);
                end
            end
            imagen_reconstruida(x,y) = real(suma) / (Nx * Ny);
        end
    end
end