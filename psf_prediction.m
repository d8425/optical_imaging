clc,clear
% make blur
img = single(imread('G:\\HD160--\\21\\G31.MIM'));

%  PSF 点
delta = 2;
xc = 6;
yc = 6;
xl = 11;
yl = 11;
psf_f = zeros([xl, yl]);
for i = 1:xl
    for j = 1:yl
        x = i;
        y = j;
        psf_f(i, j) = (1/(2*pi*(delta^2))) * exp(-((((x - xc)^2) + ((y - yc)^2))/(2*(delta^2))));
    end
end


img = conv2(img, psf_f);

%% solve psf
% blur_img = imread('');

size_range = [2,21];
sigma = [1,9];
step_size = 1;
step_sigma = 1;


epoch = 1;

reslut = [];
temp_result = 0;

[seed_size,sigma_seed] = (rand_seed(size_range,sigma));
seed_size = round(seed_size);
sigma_seed = round(sigma_seed);
for i = 1:epoch
    % if i==1
    %     size_num = seed_size;
    %     sigma_num = sigma_seed;
    % end
    % ifft_img = fft_ifft(img,sigma_num,size_num);
    % temp_result = criter(ifft_img);
    % reslut(i) = temp_result;
    % disp(num2str(temp_result))

    % rand * 10
    start_size = seed_size;
    start_sigma = sigma_seed;
    for j = 1:10
        [seed_size,sigma_seed] = (rand_seed(size_range,sigma));
        seed_size = round(seed_size);
        sigma_seed = round(sigma_seed);
        size_num(j) = seed_size;
        sigma_num(j) = sigma_seed;
        ifft_img = fft_ifft(img,sigma_num(j),size_num(j));
        temp_result = criter(ifft_img);
        reslut(i) = temp_result;
        disp(num2str(temp_result))
    end
end

function [seed_size,sigma_seed] = rand_seed(size_range,sigma)
    seed_size = rand(1)*(size_range(2)-size_range(1))+size_range(1);
    sigma_seed = rand(1)*(sigma(2)-sigma(1))+sigma(1);
end

function ifft_img = fft_ifft(img,sigma,size_n)
    xc = floor(size_n/2)+1;
    yc = floor(size_n/2)+1;
    xl = size_n;
    yl = size_n;
    psf_f = zeros([xl, yl]);
    for i = 1:xl
        for j = 1:yl
            x = i;
            y = j;
            psf_f(i, j) = (1/(2*pi*(sigma^2))) * exp(-((((x - xc)^2) + ((y - yc)^2))/(2*(sigma^2))));
        end
    end
    
    
    img = conv2(img, psf_f);
    
    % 计算填充量
    pad_rows = size(img, 1) - size(psf_f, 1);
    pad_cols = size(img, 2) - size(psf_f, 2);
    pad_top = floor(pad_rows / 2);
    pad_bottom = pad_rows - pad_top;
    pad_left = floor(pad_cols / 2);
    pad_right = pad_cols - pad_left;
    
    padded_kernel = padarray(psf_f, [pad_top, pad_left], 0, 'post');
    padded_kernel = padarray(padded_kernel, [pad_bottom, pad_right], 0, 'pre');
    
    fft_img = fft2(img);
    fft_kernel = fft2(padded_kernel);
    
    epsn = 1e-10;
    fft_deblurred = fft_img ./ (fft_kernel + epsn);
    
    deblurred_img = ifft2(fft_deblurred);
    
    deblurred_img = real(deblurred_img);
    
    deblurred_img = max(0, min(deblurred_img, 65535));
    deblurred_img = uint16(deblurred_img);
    
    ifft_img = fftshift(deblurred_img);
end

function result = criter(img_pred) % or (img_pred,img_gt) for image
    [temp, ~] =  imgradient(single(img_pred));
    result = mean(temp(:));
end
