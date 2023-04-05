%Fit instrument response function

%clears workspace
clear

path = 'C:\Users\an19697\OneDrive - University of Bristol\Desktop\Flavins\TCSPC\20220714\'
name = 'IRF_2'
path_name = [path name]
%imports dataset in current folder and specified title
data = dlmread(path_name);

%separates raw data into two matrices
t = data(2270:2550,1); %20760:20830 2400:2700
dat = data(2270:2550,2);
plot(data(2270:2550,1),data(2270:2550,2))
%converts ps to ns
t = t./1000;

%normalises data with subscript normv.m
dat = normv(dat);

[M,I]=max(dat);

%resets time zero
t0 = t(I);
t = t-t0;

%208300
%finds indices of t=0
t0ind = find(t ==0);

% isolates IRF
%t = t(t0ind-500:t0ind+500);
%dat = dat(t0ind-500:t0ind+500);

plot(t,dat);

%defines data in terms of x and y
x = t;
y = dat-0.008;

%gaussian distribution as defined by Mathematica: https://mathworld.wolfram.com/GaussianFunction.html
% Where the FWHM is given as: 2(2ln2)^0.5*sigma. Sigma is what we need to
% feed into the next step
func = @(A,t0,sigma,x) ((A./(sigma.*(2*pi).^0.5)).*exp(-((x-t0).^2)./(2.*sigma.^2)));

[curvefit,gof] = fit(x,y,func,'StartPoint',[1,0,0.14],'Robust', 'Bisquare','DiffMaxChange',0.0001,'TolFun',1.0e-15,'TolX',1e-15)

%retrieve fit coefficients
coeffs = coeffvalues(curvefit);
A = coeffs(1);
t0 = coeffs(2);
sigma = coeffs(3);

%plots fit overlaid with IRF
hold on
plot(x,func(A,t0,sigma,x),'r')
hold off

%value we need to fit data
sigma