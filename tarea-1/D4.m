function IM_GRAPPA = grappa_reconstruccion(C, m, R)
%GRAPPA_RECONSTRUCCION Reconstrucción GRAPPA para submuestreo Rx
%
% Entradas:
%   C : mapas de sensibilidad [Nx, Ny, Nc]
%   m : imagen de referencia [Nx, Ny]
%   R : factor de aceleración (ej. 2)
%
% Salida:
%   IM_GRAPPA : imagen final combinada mediante sensibilidad ponderada

    [Nx, Ny, Nc] = size(C);

    %% 1. Generar imágenes por bobina
    rec_by_coil = repmat(m, 1, 1, Nc) .* C;

    figure('Name','Imágenes originales por bobina');
    for c = 1:Nc
        subplot(2,4,c);
        imshow(abs(rec_by_coil(:,:,c)), []);
        title(['Bobina ', num2str(c)]);
    end

    %% 2. Transformada a k-space
    kspace_full = zeros(Nx, Ny, Nc);
    for c = 1:Nc
        kspace_full(:,:,c) = fftshift(fft2(rec_by_coil(:,:,c)));
    end

    %% 3. Submuestreo uniforme + líneas ACS
    nACS = 24;  % líneas centrales de calibración
    mask = zeros(Nx, Ny);

    % líneas aceleradas
    mask(1:R:end, :) = 1;

    % región ACS central
    acs_start = floor(Nx/2 - nACS/2) + 1;
    acs_end   = acs_start + nACS - 1;
    mask(acs_start:acs_end, :) = 1;

    % aplicar máscara
    kspace_us = zeros(Nx, Ny, Nc);
    for c = 1:Nc
        kspace_us(:,:,c) = kspace_full(:,:,c) .* mask;
    end

    %% 4. Mostrar k-space submuestreado
    figure('Name','K-space submuestreado');
    for c = 1:Nc
        subplot(2,4,c);
        imshow(log(abs(kspace_us(:,:,c)) + 1), []);
        title(['k-space Bobina ', num2str(c)]);
    end

    %% 5. Calibración GRAPPA (kernel 2x3)
    kernel_pe = 2;   % dirección de fase
    kernel_fe = 3;   % dirección de frecuencia

    ACS = kspace_full(acs_start:acs_end,:,:);

    src = [];
    tgt = [];

    for y = 2:size(ACS,1)-1
        for x = 2:Ny-1

            % usar solo líneas adquiridas
            if mod(y-1,R) == 0
                source_block = ACS([y-1 y+1], x-1:x+1, :);
                target_point = ACS(y, x, :);

                src = [src; source_block(:).'];
                tgt = [tgt; target_point(:).'];
            end
        end
    end

    % Pesos GRAPPA
    W = pinv(src) * tgt;

    %% 6. Interpolación de líneas faltantes
    kspace_grappa = kspace_us;

    for y = 2:Nx-1
        if mask(y,1) == 0
            for x = 2:Ny-1
                source_block = kspace_grappa([y-1 y+1], x-1:x+1, :);
                kspace_grappa(y,x,:) = reshape(source_block(:).'*W, [1 1 Nc]);
            end
        end
    end

    %% 7. Reconstrucción por bobina
    coil_recon = zeros(Nx, Ny, Nc);

    figure('Name','Reconstrucción GRAPPA por bobina');
    for c = 1:Nc
        coil_recon(:,:,c) = ifft2(ifftshift(kspace_grappa(:,:,c)));

        subplot(2,4,c);
        imshow(abs(coil_recon(:,:,c)), []);
        title(['GRAPPA Bobina ', num2str(c)]);
    end

    %% 8. Combinación ponderada por sensibilidad
    numerador = sum(conj(C) .* coil_recon, 3);
    denominador = sqrt(sum(abs(C).^2, 3));
    denominador(denominador == 0) = eps;

    IM_GRAPPA = abs(numerador ./ denominador);
    IM_GRAPPA = IM_GRAPPA / max(IM_GRAPPA(:));

    %% 9. Mostrar resultado final
    figure('Name','Reconstrucción final GRAPPA');
    subplot(1,3,1);
    imshow(m, []);
    title('Imagen original');

    subplot(1,3,2);
    imshow(IM_GRAPPA, []);
    title(['GRAPPA x', num2str(R)]);

    subplot(1,3,3);
    imshow(abs(m - IM_GRAPPA), []);
    title('Error absoluto');
end