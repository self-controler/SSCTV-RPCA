%%测试那个块放大的程序
clc;clear;close all;
X = imread('SS.png');
X = double(X)/256;
subplot(221);imshow(X);title('Orl');
Y1 = WindowGig(X,[0.45,0.2],0.2);
subplot(222);imshow(Y1);title('windowBig');
X2 = rgb2gray(X);
subplot(223);imshow(X2);title('Orl');
Y2 = WindowGig(X2,[0.45,0.2],0.2);
subplot(224);imshow(Y2);title('windowBig');
% imwrite(X2,'SSgry.png');
% imwrite(Y1,'windBig.png');
% imwrite(Y2,'gryWindBig.png');