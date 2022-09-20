clc;
clear;
close all;

%% 读取图像序列
path = 'images\seq1\';
imgDir = dir([path '*.bmp']);
len = length(imgDir);

for i = 1:len
    temp = imread([path, num2str(i),'.bmp']);
    img(:,:,i) = temp;
end
img = double(img);
[m,n,l] = size(img);

%% 滤波器设计
fs=50;
% 十帧累加滤波器
frame = 10;
b = cat(2,1,zeros(1,frame-2),-1);
a = cat(2,frame,-frame,zeros(1,frame-2));
[h,f] = freqz(b,a,1000,fs);
figure;plot(f,20*log10(abs(h)),'linewidth',2);hold on;
xlabel('Normalized Frequency (f/fs)');ylabel('Magnitude (dB)');

% 0.1hz陷波滤波器
b01 = [1 -2*cos(2*pi*0.1/fs) 1];
a01 = [1 -2*cos(2*pi*0.1/fs) 0.999999^2];
[h,f] = freqz(b01,a01,1000,fs);
plot(f,20*log10(abs(h)),'linewidth',2);hold on;

% 1hz陷波滤波器
b1 = [1 -2*cos(2*pi*1/fs) 1];
a1 = [1 -2*cos(2*pi*1/fs) 0.999^2];
[h,f] = freqz(b1,a1,1000,fs);
plot(f,20*log10(abs(h)),'linewidth',2);hold on;

% 5.5hz陷波滤波器
b55 = [1 -2*cos(2*pi*5.5/fs) 1];
a55 = [1 -2*cos(2*pi*5.5/fs) 0.999^2];
[h,f] = freqz(b55,a55,1000,fs);
plot(f,20*log10(abs(h)),'linewidth',2);hold on;

% 11.1hz陷波滤波器
b111 = [1 -2*cos(2*pi*11.1/fs) 1];
a111 = [1 -2*cos(2*pi*11.1/fs) 0.999^2];
[h,f] = freqz(b111,a111,1000,fs);hold on;
plot(f,20*log10(abs(h)),'linewidth',2);

legend('多帧累加滤波器','0.1hz陷波','1hz陷波','5.5hz陷波','11.1hz陷波')

%% 10帧累加
for i = 1:m
    for j = 1:n
        img_accum(i,j,:)=filter(b,a,img(i,j,:));
    end
end
diff = img_accum(:,:,end)-img(:,:,end);
diff = mat2gray(diff);

figure;
subplot(131);imshow(uint8(img(:,:,end)));title('原图')
subplot(132);imshow(uint8(img_accum(:,:,end)));title('滤波器版的10帧累加')
subplot(133);imshow(diff);title('差分')

%% 对每个像素点做0.1hz陷波
for i = 1:m
    for j = 1:n
        img_01hz(i,j,:)=filter(b01,a01,img(i,j,:));
    end
end
diff = img_01hz(:,:,end)-img(:,:,end);
diff = mat2gray(diff);

figure;
subplot(131);imshow(uint8(img(:,:,end)));title('原图')
subplot(132);imshow(uint8(img_01hz(:,:,end)));title('0.1hz陷波后')
subplot(133);imshow(diff);title('差分')

%% 对每个像素点做1hz陷波
for i = 1:m
    for j = 1:n
        img_1hz(i,j,:)=filter(b1,a1,img(i,j,:));
    end
end
diff = img_1hz(:,:,end)-img(:,:,end);
diff = mat2gray(diff);

figure;
subplot(131);imshow(uint8(img(:,:,end)));title('原图')
subplot(132);imshow(uint8(img_1hz(:,:,end)));title('1hz陷波后')
subplot(133);imshow(diff);title('差分')

%% 对每个像素点做5.5hz陷波
for i = 1:m
    for j = 1:n
        img_55hz(i,j,:)=filter(b55,a55,img(i,j,:));
    end
end
diff = img_55hz(:,:,end)-img(:,:,end);
diff = mat2gray(diff);

figure;
subplot(131);imshow(uint8(img(:,:,end)));title('原图')
subplot(132);imshow(uint8(img_55hz(:,:,end)));title('5.5hz陷波后')
subplot(133);imshow(diff);title('差分')

%% 对每个像素点做11.1hz陷波
for i = 1:m
    for j = 1:n
        img_111hz(i,j,:)=filter(b111,a111,img(i,j,:));
    end
end
diff = img_111hz(:,:,end)-img(:,:,end);
diff = mat2gray(diff);

figure;
subplot(131);imshow(uint8(img(:,:,end)));title('原图')
subplot(132);imshow(uint8(img_111hz(:,:,end)));title('11.1hz陷波后')
subplot(133);imshow(diff);title('差分')
