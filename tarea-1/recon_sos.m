% recon_sos.m
% Carga `m.mat` y `C.mat`, genera la imagen por cada bobina
% y calcula la reconstrucción por suma de cuadrados (SOS).

% Cargar archivos (asume que están en el mismo directorio)
data = load('m.mat');
Cdata = load('C.mat');

% Extraer variables (no asumimos nombres fijos)
fn = fieldnames(data);
m = data.(fn{1});
fnc = fieldnames(Cdata);
C = Cdata.(fnc{1});

% Asegurar dimensiones: `m` debe ser 2D, `C` debe ser Nx x Ny x Nc
if ndims(m) ~= 2
    error('La variable m no tiene dimensiones 2D esperadas');
end

% Si C viene como 2D con cada bobina en columnas, intentar reconstruir
if ndims(C) == 2 && size(C,1) == numel(m)
    % C es (Nx*Ny) x Nc
    [Nx,Ny] = size(m);
    Nc = size(C,2);
    C = reshape(C, [Nx, Ny, Nc]);
end

[Nx,Ny,maybeNc] = size(C);
if maybeNc == 1 && isequal([Nx,Ny], size(m))
    Nc = 1;
elseif isequal([Nx,Ny], size(m))
    Nc = maybeNc;
else
    % Si no coincide, intentar invertir dimensiones
    [Nx2,Ny2,Nc2] = size(C);
    if Nx2==size(m,2) && Ny2==size(m,1)
        C = permute(C, [2 1 3]);
        Nc = Nc2;
    else
        error('Dimensiones de C no coinciden con m');
    end
end

% Inicializar matriz para imágenes por bobina
coil_images = zeros(size(m,1), size(m,2), Nc);

for k = 1:Nc
    coil_images(:,:,k) = m .* C(:,:,k);
end

% Reconstrucción por suma de cuadrados
sos = sqrt(sum(abs(coil_images).^2, 3));

% Mostrar resultados
figure('Name','Imágenes por bobina');
nplot = min(Nc, 16);
cols = ceil(sqrt(nplot));
rows = ceil(nplot/cols);
for k = 1:nplot
    subplot(rows, cols, k)
    imagesc(abs(coil_images(:,:,k)))
    colormap gray
    axis image off
    title(sprintf('Bobina %d', k))
end

figure('Name','Reconstrucción SOS');
imagesc(abs(sos)); colormap gray; axis image off
title('Reconstrucción por Suma de Cuadrados (SOS)');

% Guardar resultados
save('results_sos.mat', 'coil_images', 'sos');

fprintf('Hecho: guardado en results_sos.mat (coil_images, sos)\n');
