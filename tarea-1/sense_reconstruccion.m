function IM_R = sense_reconstruccion(C, m, rate)
%SENSE_RECONSTRUCCION Reconstrucción SENSE para submuestreo Rx
%
% Entradas:
%   C    : mapas de sensibilidad [Nx, Ny, Nc]
%   m    : imagen original [Nx, Ny]
%   rate : factor de aceleración (2, 4, 6, 8, ...)
%
% Salida:
%   IM_R : imagen reconstruida y normalizada

    % Dimensiones
    [Nx, Ny, Nc] = size(C);

    % Generar imágenes por bobina
    rec_by_coil = repmat(m, 1, 1, Nc) .* C;

    % Mostrar imágenes por bobina
    figure('Name', ['Imágenes por bobina - x', num2str(rate)]);
    for n = 1:Nc
        subplot(2,4,n);
        imshow(abs(rec_by_coil(:,:,n)), []);
        title(['Bobina ', num2str(n)]);
    end

    % Máscara de submuestreo
    mask = zeros(size(C));
    mask(1:rate:end,:,:) = 1;

    % K-space submuestreado
    kspace = zeros(size(C));
    figure('Name', ['K-space submuestreado x', num2str(rate)]);
    for n = 1:Nc
        kspace(:,:,n) = fftshift(fft2(rec_by_coil(:,:,n))) .* mask(:,:,n);
        subplot(2,4,n);
        imshow(log(abs(kspace(:,:,n)) + 1), []);
        title(['k-space bobina ', num2str(n)]);
    end

    % Extraer solo líneas adquiridas
    idx_sampled = 1:rate:Nx;
    M = numel(idx_sampled);

    brain_sub = complex(zeros(M, Ny, Nc));
    kspace_sub = complex(zeros(M, Ny, Nc));

    figure('Name', ['Bobinas submuestreadas x', num2str(rate)]);
    for n = 1:Nc
        kspace_sub(:,:,n) = kspace(idx_sampled,:,n);
        brain_sub(:,:,n) = ifft2(fftshift(kspace_sub(:,:,n)));

        subplot(2,4,n);
        imshow(abs(brain_sub(:,:,n)), []);
        title(['Bobina ', num2str(n), ' x', num2str(rate)]);
    end

    % Reconstrucción SENSE
    IMSENSE = complex(zeros(Nx, Ny));

    for jj = 1:M
        idx_folded = jj:rate:Nx;

        for ii = 1:Ny
            % Matriz de sensibilidades
            CC = squeeze(C(idx_folded, ii, :)).';

            % Señales medidas en las bobinas
            a = squeeze(brain_sub(jj, ii, :));

            % Solución por pseudo-inversa
            IMSENSE(idx_folded, ii) = pinv(CC) * a;
        end
    end

    % Magnitud y normalización
    IMSENSE = abs(IMSENSE);
    IM_R = IMSENSE / max(IMSENSE(:));

    % Mostrar resultado final
    figure('Name', ['Reconstrucción SENSE x', num2str(rate)]);
    subplot(1,3,1);
    imshow(m, []);
    title('Imagen original');

    subplot(1,3,2);
    imshow(IM_R, []);
    title(['Reconstrucción SENSE x', num2str(rate)]);

    subplot(1,3,3);
    imshow(abs(m - IM_R), []);
    title('Error absoluto');
end