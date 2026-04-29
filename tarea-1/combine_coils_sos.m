function img_sos = combine_coils_sos(coil_imgs)
%COMBINE_COILS_SOS Combina imágenes multibobina usando
% el método Sum of Squares (SoS).
%
% Entrada:
%   coil_imgs : arreglo [Nx, Ny, Nc]
%               Imágenes complejas de cada bobina.
%
% Salida:
%   img_sos   : imagen reconstruida [Nx, Ny]

    % Suma de cuadrados de la magnitud
    img_sos = sqrt(sum(abs(coil_imgs).^2, 3));
end