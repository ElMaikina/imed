% RECONSTRUCCION_SENSE
% cargamos la imagen y los mapas de sensibilidad de las bobinas
load('IMAGENES_CURSO\DATOS_SENSE\C.mat')
load('IMAGENES_CURSO\DATOS_SENSE\m.mat')

figure, 
imshow(m,[])
title('Imagen referencia')

% multiplicamos la imagen original por cada uno de los mapas de
% sensibilidad de las bobinas
rec_by_coil = repmat(m,1,1,size(C,3)).*C;

figure, 
for n=1:size(C,3)
    subplot(2,4,n), imshow(abs(rec_by_coil(:,:,n)),[]), title(['bobina ',num2str(n)])
end

% generamos un submuestreo uniforme en la direccion de las fases
mask = zeros(size(C));
mask(1:2:end,:,:) = mask(1:2:end,:,:) + 1;
kspace = [];

figure,
for n=1:size(C,3)
    kspace(:,:,n) = fftshift(fft2(rec_by_coil(:,:,n))).*mask(:,:,n);
    subplot(2,4,n), imshow(log(abs(kspace(:,:,n))),[]), title(['kspace bobina x2',num2str(n)])
end

% mostramos la imagen resultante de cada una de las bobinas
im_rec_fullsize = [];
figure,
for n=1:size(C,3)
    im_rec_fullsize(:,:,n) = abs(ifft2(fftshift(kspace(:,:,n))));
    subplot(2,4,n), imshow(im_rec_fullsize(:,:,n),[]), title(['im-recon bobina x2',num2str(n)])
end

% mostramos la imagen resultante de cada una de las bobinas considerando
% solo las lineas con informacion. 
kspace_subx2 = [];
im_rec_fullsize = [];
figure,
for n=1:size(C,3)
    kspace_subx2(:,:,n) = kspace(1:2:end,:,n);
    brain_subx2(:,:,n) = ifft2(fftshift(kspace_subx2(:,:,n)));
    subplot(4,4,n), imshow(abs(brain_subx2(:,:,n)),[]), title(['im-recon-sub bobina x2',num2str(n)])
    subplot(4,4,n+8), imshow(log(abs(kspace_subx2(:,:,n))),[]), title(['kspace-sub bobina x2',num2str(n)])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% para sense con submuestreo x2
M = size(brain_subx2,1);
N = size(brain_subx2,2);
Ncoils = 8 ;
rate = 2; %submuestreo
MM = size(C,1);

% para X2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% #################### SENSE ####################################################
IMSENSE = complex(zeros([MM,N]));
for jj=1:M % rows
    for ii=1:N % columns
        % considerar que la traspuesta compleja es ".'"
        CC = squeeze(C([jj , jj + M], ii , : )).';
        a = squeeze(brain_subx2(jj,ii,:));
        IMSENSE([jj , jj + M] , ii ) = inv(CC.'*inv(eye(size(CC,1)))*CC)*CC.'*inv(eye(size(CC,1)))*a;
    end
end

IMSENSE = abs(IMSENSE);
IM_R = IMSENSE/max(IMSENSE(:));

figure,
subplot 131, imshow(m)
title('brain-original')                                       
subplot 132, imshow(IM_R)
title('brain-recon b0x2')                                    
subplot 133, imshow(abs(m-IM_R))
title('diff')                                    