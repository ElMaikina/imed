% Funcion que interpola una imagen 3D a una salida de diferentes
% dimensiones a partir de un factor dado. La implementacion esta
% hecha a partir del metodo indicado en:
%
%   - https://es.wikipedia.org/wiki/Interpolaci%C3%B3n_trilineal#cite_note-NASA-4
%
function volumen_salida = trilineal_MS_SS(volumen_entrada, factor)
    % Se revisa que el factor cumpla con el enunciado
    if factor <= 0.05 || factor >= 2
        error('El factor debe estar en el rango (0.05, 2)');
    end
    % Obtiene las dimensiones de la matriz de entrada para leer
    [Nx, Ny, Nz] = size(volumen_entrada);
    % Genera el espacio de salida para escribir
    Nx_out = round(Nx * factor);
    Ny_out = round(Ny * factor);
    Nz_out = round(Nz * factor);
    % Rellena la matriz de salida con ceros
    volumen_salida = zeros(Nx_out, Ny_out, Nz_out);
    % Convierte las entradas a doble precision
    volumen_entrada = double(volumen_entrada);
    % Itera a partir de las tres dimensiones para mapear los puntos
    for x_out = 1:Nx_out
        for y_out = 1:Ny_out
            for z_out = 1:Nz_out
                % Lee un punto en la matriz de entrada
                x_in = (x_out - 1) / factor + 1;
                y_in = (y_out - 1) / factor + 1;
                z_in = (z_out - 1) / factor + 1;
                % Obtiene los puntos vecinos anteriores aproximando "hacia abajo"
                x0 = floor(x_in);
                y0 = floor(y_in);
                z0 = floor(z_in);
                % Obtiene los puntos vecinos posteriores aproximando "hacia arriba"
                x1 = x0 + 1;
                y1 = y0 + 1;
                z1 = z0 + 1;
                % Usa clamping para manejar los bordes de la matriz
                x0 = max(1, min(x0, Nx));
                x1 = max(1, min(x1, Nx));
                y0 = max(1, min(y0, Ny));
                y1 = max(1, min(y1, Ny));
                z0 = max(1, min(z0, Nz));
                z1 = max(1, min(z1, Nz));
                % Calcula las distancias relativas (pesos)
                dx = x_in - x0;
                dy = y_in - y0;
                dz = z_in - z0;
                % Obtiene los ocho vecinos que rodean nuestra coordenada central
                I000 = volumen_entrada(x0, y0, z0);
                I100 = volumen_entrada(x1, y0, z0);
                I010 = volumen_entrada(x0, y1, z0);
                I110 = volumen_entrada(x1, y1, z0);
                I001 = volumen_entrada(x0, y0, z1);
                I101 = volumen_entrada(x1, y0, z1);
                I011 = volumen_entrada(x0, y1, z1);
                I111 = volumen_entrada(x1, y1, z1);
                % Interpola en el eje "X"
                c00 = (1 - dx) * I000 + dx * I100;
                c01 = (1 - dx) * I001 + dx * I101;
                c10 = (1 - dx) * I010 + dx * I110;
                c11 = (1 - dx) * I011 + dx * I111;
                % Interpola en el eje "Y"
                c0 = (1 - dy) * c00 + dy * c10;
                c1 = (1 - dy) * c01 + dy * c11;
                % Interpola en el eje "Z"
                valor_interpolado = (1 - dz) * c0 + dz * c1;
                % Finalmente guarda el punto en la matriz de salida
                volumen_salida(x_out, y_out, z_out) = valor_interpolado;

            end
        end
    end
end