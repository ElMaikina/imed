function imagen_reconstruida = reconstruccion_MS_SS(imagen)
    % Convertir a double para precisión
    imagen = double(imagen);
    % Obtener dimensiones
    [Nx, Ny] = size(imagen);
    %--------------------------------------------------
    % 1. Calcular DFT 2D manual
    %--------------------------------------------------
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
    %--------------------------------------------------
    % 2. Crear máscara de filtrado
    %--------------------------------------------------
    mascara = ones(Nx, Ny);
    
    centro_x = round(Nx/2);
    centro_y = round(Ny/2);
    
    %----------------------------------------------
    % Caso 1: eliminar puntos periódicos
    %----------------------------------------------
    % Eliminamos puntos equiespaciados en el espectro
    
    separacion = 10; % ajustar según imagen
    
    for u = 1:separacion:Nx
        for v = 1:separacion:Ny
            mascara(u,v) = 0;
        end
    end
    
    %----------------------------------------------
    % Caso 2: eliminar línea diagonal (franjas)
    %----------------------------------------------
    for u = 1:Nx
        for v = 1:Ny
            % condición de diagonal
            if abs((u - centro_x) - (v - centro_y)) < 3
                mascara(u,v) = 0;
            end
        end
    end
    
    %--------------------------------------------------
    % 3. Aplicar filtro
    %--------------------------------------------------
    espectro_filtrado = espectro .* mascara;
    
    %--------------------------------------------------
    % 4. IDFT 2D manual
    %--------------------------------------------------
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