f = -4:4;
t = -4:4; % Zaman aralığı
x = zeros(size(t)); % x(t) sinyali için boş bir dizi oluşturuyoruz
% -4<t<-2 aralığındaki değerler
x(t > -4 & t < -2) = -1;
% -2<t<2 aralığındaki değerler
x(t > -2 & t < 2) = 0.5 * t(t > -2 & t < 2);
% 2<t<4 aralığındaki değerler
x(t > 2 & t < 4) = 1;

% Fourier dönüşümünü hesaplıyoruz
Xf = fft(x);

% Fourier dönüşümünü çizdiriyoruz
figure;
plot(f, abs(fftshift(Xf)));  %Fourier dönüşüm grafiğini çizdiriyoruz.
title('x(t) Fourier Dönüşümü');
xlabel('Frekans (Hz)');
ylabel('|X(f)|');

figure;
subplot(3,1,1); %Fonksiyonların grafiğini görebilmek için çıktı ekranını 3e böler ve 1. konuma yerleştirir.
plot(f,real(Xf),f,imag(Xf)); % Xf fonksiyonunun grafiği çizdirilir.
legend('Real Part', 'Imaginary Part'); % Oluşan grafikler isimlendirilir.

%x fonksiyonu için x ve y eksenlerine isim atamaları yapılır ve grafiğe başlık eklenir.
xlabel('f'); 
ylabel('X(f) function');
title('X(f)');

subplot(3,1,2); %Fonksiyonların grafiğini görebilmek için çıktı ekranını 3e böler ve 2. konuma yerleştirir.
plot(f,abs(Xf)); % genlik grafiği için mutlak fonksiyonu kullanılır.
xlabel('f');
ylabel('Genlik');
title('X(f) Genlik Grafiği');

subplot(3,1,3); %Fonksiyonların grafiğini görebilmek için çıktı ekranını 3e böler ve 3. konuma yerleştirir.
plot(f,angle(Xf)); % faz grafiği için açı fonksiyonu kullanılır.
xlabel('f');
ylabel('faz');
title('X(f) faz Grafiği');



% e^(-j2pift0), t0: (t-t0) öteleme miktarı.
fb = 1:9;
tb = 1:9; % Zaman aralığı
xb = zeros(size(tb)); % x(t) sinyali için boş bir dizi oluşturuyoruz

% 1<t<3 aralığındaki değerler
xb(tb > 1 & tb < 3) = -1;
% 3<t<7 aralığındaki değerler
xb(tb > 3 & tb < 7) = 0.5 * tb(tb > 3 & tb < 7);
% 7<t<9 aralığındaki değerler
xb(tb > 7 & tb < 9) = 1;

% Fourier dönüşümünü hesaplıyoruz
Xfb = fft(xb);
figure; %Yeni bir figure ekranı açtırıyoruz.
plot(fb, abs(fftshift(Xfb))); %Fourier dönüşüm grafiğini çizdiriyoruz.
title('x(t-5) Fourier Dönüşümü');
xlabel('Frekans (Hz)');
ylabel('|X(f)|');

figure;
subplot(3,1,1); %Fonksiyonların grafiğini görebilmek için çıktı ekranını 3e böler ve 1. konuma yerleştirir.
plot(fb,real(Xfb),fb,imag(Xfb)); % Xf fonksiyonunun grafiği çizdirilir.
legend('Real Part', 'Imaginary Part'); % Oluşan grafikler isimlendirilir.

%x fonksiyonu için x ve y eksenlerine isim atamaları yapılır ve grafiğe başlık eklenir.
xlabel('f'); 
ylabel('X(fb) function');
title('X(fb)');

subplot(3,1,2); %Fonksiyonların grafiğini görebilmek için çıktı ekranını 3e böler ve 2. konuma yerleştirir.
plot(fb,abs(Xfb)); % genlik grafiği için mutlak fonksiyonu kullanılır.
xlabel('f');
ylabel('Genlik');
title('X(fb) Genlik Grafiği');

subplot(3,1,3); %Fonksiyonların grafiğini görebilmek için çıktı ekranını 3e böler ve 3. konuma yerleştirir.
plot(fb,angle(Xfb)); % faz grafiği için açı fonksiyonu kullanılır.
xlabel('f');
ylabel('faz');
title('X(fb) faz Grafiği');



x2 = exp(1i * 10 * pi * t); % İşaretin tanımı
Xf2 = fft(x2); %İşaretin fourier dönüşümü
% f,t,x ve Xf en üstte tanımlı olduğu ve aynı değerlerden oluşacağı için tekrardan tanımlamıyoruz. 

% Fourier dönüşümünü konvolüsyonla hesaplıyoruz
Xfc= conv(Xf2,Xf, 'same') * 0.01; % 0.01 ile örnekleme periyodu çarpılır

% Fourier dönüşümünü çizdiriyoruz
figure;
plot(f, abs(fftshift(Xfc)));  %Fourier dönüşüm grafiğini çizdiriyoruz.
title('xc(t) Fourier Dönüşümü');
xlabel('Frekans (Hz)');
ylabel('|Xc(f)|');

figure;
subplot(3,1,1); %Fonksiyonların grafiğini görebilmek için çıktı ekranını 3e böler ve 1. konuma yerleştirir.
plot(f,real(Xfc),f,imag(Xfc)); % Xf fonksiyonunun grafiği çizdirilir.
legend('Real Part', 'Imaginary Part'); % Oluşan grafikler isimlendirilir.

%x fonksiyonu için x ve y eksenlerine isim atamaları yapılır ve grafiğe başlık eklenir.
xlabel('f'); 
ylabel('Xc(f) function');
title('Xc(f)');

subplot(3,1,2); %Fonksiyonların grafiğini görebilmek için çıktı ekranını 3e böler ve 2. konuma yerleştirir.
plot(f,abs(Xfc)); % genlik grafiği için mutlak fonksiyonu kullanılır.
xlabel('f');
ylabel('Genlik');
title('Xc(f) Genlik Grafiği');

subplot(3,1,3); %Fonksiyonların grafiğini görebilmek için çıktı ekranını 3e böler ve 3. konuma yerleştirir.
plot(f,angle(Xfc)); % faz grafiği için açı fonksiyonu kullanılır.
xlabel('f');
ylabel('faz');
title('Xc(f) faz Grafiği');


f3 = -(4/3):(4/3);
t3 = -(4/3):(4/3); % Zaman aralığı ölçeklemeyi yapıyoruz
xd = zeros(size(t3)); % x(t) sinyali için boş bir dizi oluşturuyoruz
% -4<t<-2 aralığındaki değerler
xd(t3 > -4/3 & t3 < -2/3) = -1;
% -2<t<2 aralığındaki değerler
xd(t3 > -2/3 & t3 < 2/3) = 0.5 * t3(t3 > -2/3 & t3 < 2/3);
% 2<t<4 aralığındaki değerler
xd(t3 > 2/3 & t3 < 4/3) = 1;
% x(t) fonksiyonunun ölçeklendirilmiş halini hesaplayalım.
scaled_x = (1/3) * xd;

Xfscaled = fft(scaled_x); %İşaretin fourier dönüşümü
% Fourier dönüşümünü çizdiriyoruz
figure;
plot(f3, abs(fftshift(Xfscaled)));  %Fourier dönüşüm grafiğini çizdiriyoruz.
title('xscaled(t) Fourier Dönüşümü');
xlabel('Frekans (Hz)');
ylabel('|Xscaled(f)|');

figure;
subplot(3,1,1); %Fonksiyonların grafiğini görebilmek için çıktı ekranını 3e böler ve 1. konuma yerleştirir.
plot(f3,real(Xfscaled),f3,imag(Xfscaled)); % Xf fonksiyonunun grafiği çizdirilir.
legend('Real Part', 'Imaginary Part'); % Oluşan grafikler isimlendirilir.

%x fonksiyonu için x ve y eksenlerine isim atamaları yapılır ve grafiğe başlık eklenir.
xlabel('f'); 
ylabel('Xscaled(f) function');
title('Xscaled(f)');

subplot(3,1,2); %Fonksiyonların grafiğini görebilmek için çıktı ekranını 3e böler ve 2. konuma yerleştirir.
plot(f3,abs(Xfscaled)); % genlik grafiği için mutlak fonksiyonu kullanılır.
xlabel('f');
ylabel('Genlik');
title('Xscaled(f) Genlik Grafiği');

subplot(3,1,3); %Fonksiyonların grafiğini görebilmek için çıktı ekranını 3e böler ve 3. konuma yerleştirir.
plot(f3,angle(Xfscaled)); % faz grafiği için açı fonksiyonu kullanılır.
xlabel('f');
ylabel('faz');
title('Xscaled(f) faz Grafiği');


% 𝑦(𝑡) = 𝑠𝑔𝑛(𝑡)
% f ve t en üstte tanımlı olduğu ve aynı değerlerden oluşacağı için tekrardan tanımlamıyoruz.
xsgn = sign(t);
% Fourier dönüşümünü hesaplıyoruz
Xfsgn = fft(xsgn);

% Fourier dönüşümünü çizdiriyoruz
figure;
plot(f, abs(fftshift(Xfsgn)));  %Fourier dönüşüm grafiğini çizdiriyoruz.
title('x(t) Fourier Dönüşümü');
xlabel('Frekans (Hz)');
ylabel('|Xsgn(f)|');

figure;
subplot(3,1,1); %Fonksiyonların grafiğini görebilmek için çıktı ekranını 3e böler ve 1. konuma yerleştirir.
plot(f,real(Xfsgn),f,imag(Xfsgn)); % Xf fonksiyonunun grafiği çizdirilir.
legend('Real Part', 'Imaginary Part'); % Oluşan grafikler isimlendirilir.

%x fonksiyonu için x ve y eksenlerine isim atamaları yapılır ve grafiğe başlık eklenir.
xlabel('f'); 
ylabel('Xsgn(f) function');
title('Xsgn(f)');

subplot(3,1,2); %Fonksiyonların grafiğini görebilmek için çıktı ekranını 3e böler ve 2. konuma yerleştirir.
plot(f,abs(Xfsgn)); % genlik grafiği için mutlak fonksiyonu kullanılır.
xlabel('f');
ylabel('Genlik');
title('Xsgn(f) Genlik Grafiği');

subplot(3,1,3); %Fonksiyonların grafiğini görebilmek için çıktı ekranını 3e böler ve 3. konuma yerleştirir.
plot(f,angle(Xfsgn)); % faz grafiği için açı fonksiyonu kullanılır.
xlabel('f');
ylabel('faz');
title('Xsgn(f) faz Grafiği');
