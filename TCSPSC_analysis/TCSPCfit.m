%clears workspace
%clear

%enter pulse duration from fitting IRF trace
pulsewidth = 0.0671

path = ('C:\Users\an19697\OneDrive - University of Bristol\Desktop\Flavins\Paulina\July\')
name = ('design4_F+H_500LP_600SP_POL_445nm_3')

%%
%%%%PATH FOR SEQUENTIALLY ACQUIRED DATA %%%%%%%%%%%%%%%%
% 
% tmp = dlmread([path name]); 
% tmp = tmp';
% t = tmp(:,1);
% 
% 
% n = size(tmp,1);
% m = size(tmp,2)./2;
% 
% tmp2 = zeros(m,n);
% 
% for i = 1:m
% tmp3 = tmp(:,2*i);
% tmp2(i,:) = tmp3' + tmp2(i,:);
% end
% 
% dat = tmp2;
% 
% t = t(20701:22801); %22001
% 
% dat = (mean(dat,1)');
% 
% dat = dat(20701:22801);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
path_name = [path name]

%imports dataset in current folder and specified title
data = dlmread(path_name);

%plot(data(:,1),data(:,2))

% %separates raw data into two matrices

t = data(2000:7000,1);
dat = data(2000:7000,2);

%plot(t,dat)

%converts ps to ns
t = t./1000;

%normalises data with subscript normv.m
dat = normv(dat);
%plot(t,dat)

%finds index of maximum data value (I);
[M,I]=max(dat);

%resets time zero
t0 = t(I);
t = t-t0;

%finds indices of t=0
t0ind = find(t ==0);

%throws away > 10 ns pre-time zero
% t = t(t0ind-5:end);
% dat = dat(t0ind-5:end);

%plot data on linear and log y-axes

subplot(2,1,1)
plot(t,dat,'Linewidth',1, 'DisplayName', 'Design 5');set(gca,'FontSize',18,'LineWidth',1);xlabel('time / ns');ylabel('Normalised intensity');axis([-1 15 0 1.1])
subplot(2,1,2)
plot(t,dat,'Linewidth',1, 'DisplayName', 'Design 5');set(gca, 'YScale', 'log');set(gca,'FontSize',18,'LineWidth',1);xlabel('time / ns');ylabel('Normalised log(intensity)');axis([-1 15 0 1.1])


% command prompt to ask for number of exponential to fit to
prompt = {'Number of exponentials to fit: 1, 2 or 3?'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'1'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
noexp = str2double(answer{1});

%subsequent command prompt to input time constant starting guesses
if noexp == 2
    prompt = {'Starting guess of time constant 1 / ns','Starting guess of time constant 2 / ns'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'50','1000'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
guessinput(1) = str2double(answer{1});
guessinput(2) = str2double(answer{2});

else if noexp == 1
       prompt = {'Starting guess of time constant / ns'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'1000'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
guessinput = str2double(answer{1});
    else
          prompt = {'Starting guess of time constant 1 / ns','Starting guess of time constant 2 / ns', 'Starting guess of time constant 3 / ns'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'1','2','3'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
guessinput(1) = str2double(answer{1});
guessinput(2) = str2double(answer{2});
guessinput(3) = str2double(answer{3});
        
    end
end

%define data to fit in terms of x and y
x = t;
y = dat;

%fit to biexponential convolved with Gaussian IRF

if noexp == 2
    
    %creates starting guess matrix with amplitudes for both exponentials
    %equal and pulse width from prior fitting
    startguess = zeros(1,6);
    startguess(1) = 0.5;
    startguess(2) = 0.5;
    startguess(3) = pulsewidth;
    startguess(4) = 0;
    %converts time constants into rate constants consistent with convolved
    %function
    startguess(5) = 1./guessinput(1);
    startguess(6) = 1./guessinput(2);
    
    %upper bounds only constraining IRF width
    upperbounds = ([Inf,Inf,pulsewidth+0.01,Inf,Inf,Inf]);
    
    %lower bounds only constraining IRF width
    lowerbounds = ([0,0,pulsewidth-0.01,-Inf,-Inf,-Inf]);
    
    %convolved function
    func = @(B1,B2,sigma,t0,k1,k2,x) ((B1.*exp(0.5.*k1.*(k1.*sigma.^2 + 2.*t0 - 2.*x)).*(1-sqrt(1./sigma.^2).*sigma.*erf((k1.*sigma.^2 + t0 - x)./(sqrt(2).*sigma))) + B2.*exp(0.5.*k2.*(k2.*sigma.^2 + 2.*t0 - 2.*x)).*(1 - sqrt(1./sigma.^2).*sigma.*erf((k2.*sigma.^2 + t0 - x)./(sqrt(2).*sigma)))))./(2.*sqrt(1/sigma.^2).*sigma);
    
    %curve fit
    
    %tf = excludedata(x,y,'domain',[0.21 0.68]); excludes data OUTSIDE the
    %specified range
    
    [curvefit,gof] = fit(x,y,func,'upper',upperbounds, 'lower',lowerbounds, 'StartPoint',startguess, 'Robust', 'Bisquare','DiffMaxChange',0.0001,'TolFun',1.0e-15,'TolX',1e-15,'Exclude', x>0.3 & x<0.7)

    %retrieve fit coefficients
    coeffs = coeffvalues(curvefit);
    B1 = coeffs(1);
    B2 = coeffs(2);
    sigma = coeffs(3);
    t0 = coeffs(4);
    k1 = coeffs(5);
    k2 = coeffs(6);

    p = confint(curvefit);
    p = p(2,:);
    n = coeffs-p;
    n(5) = 1./(coeffs(5))-1./p(5);
    n(6) = 1./(coeffs(6))-1./p(6);

    coeffs(5) = 1./coeffs(5);
    coeffs(6) = 1./coeffs(6);

    %values
    vals = {'B1','B2','sigma','t0','tau1','tau2'};

%results table
results = [coeffs' n']



subplot(2,1,1)
hold on
plot(x,func(B1,B2,sigma,t0,k1,k2,x),'r','DisplayName','No Binding Site');axis([-1 12 0 1.1])
subplot(2,1,2)
hold on
plot(x,func(B1,B2,sigma,t0,k1,k2,x),'r','DisplayName','No Binding Site');axis([-1 12 1e-5 1.1])

    
else if noexp == 1;

    startguess = zeros(1,4);
    startguess(1) = 1;
    startguess(2) = pulsewidth;
    startguess(3) = 0;
    startguess(4) = 1./guessinput(1);
    
    %upper bounds only constraining IRF width
    upperbounds = ([Inf,pulsewidth+0.01,Inf,Inf]);
    
    %lower bounds only constraining IRF width
    lowerbounds = ([0,pulsewidth-0.01,-Inf,-Inf]);
    
    func = @(B1,sigma,t0,k1,x) (B1.*exp(0.5.*k1.*(k1.*sigma.^2 + 2.*t0 - 2.*x)).*(1 - sqrt(1./sigma.^2).*sigma.*erf((k1.*sigma.^2 + t0 - x)./(sqrt(2).*sigma))))./(2.*sqrt(1./sigma.^2).*sigma);

[curvefit,gof] = fit(x,y,func, 'upper', upperbounds, 'lower', lowerbounds, 'StartPoint', startguess,  'Robust', 'LAR','MaxFunEvals',500000,'TolX',10E-8)

%get coefficients out of the fit
coeffs = coeffvalues(curvefit);
B1 = coeffs(1);
sigma = coeffs(2);
t0 = coeffs(3);
k1 = coeffs(4);

%retrieve lower bound of confidence inverval, and then subtract from error
%to get +/- error
p = confint(curvefit);
p = p(2,:);
n = coeffs-p;

n(4) = 1./(coeffs(4))-1./p(4);

coeffs(4) = 1./coeffs(4);

%values
vals = {'B1','sigma','t0','tau1'};

%results table
results = [coeffs' n'];
%note that the decay constant here is a rate, and you want tau, so take the
%reciprocal.

subplot(2,1,1)
hold on
plot(x,func(B1,sigma,t0,k1,x),'r')
subplot(2,1,2)
hold on
plot(x,func(B1,sigma,t0,k1,x))

    else
        %creates starting guess matrix with amplitudes for both exponentials
    %equal and pulse width from prior fitting
    startguess = zeros(1,8);
    startguess(1) = 0.33;
    startguess(2) = 0.33;
    startguess(3) = 0.33;
    startguess(4) = pulsewidth;
    startguess(5) = 0;
    %converts time constants into rate constants consistent with convolved
    %function
    startguess(6) = 1./guessinput(1);
    startguess(7) = 1./guessinput(2);
    startguess(8) = 1./guessinput(3);
    
    %upper bounds only constraining IRF width
    upperbounds = ([Inf,Inf,Inf, pulsewidth+0.01,Inf,Inf, Inf, Inf]);
    
    %lower bounds only constraining IRF width
    lowerbounds = ([0,0,0,pulsewidth-0.01,-Inf,0,0,0]);
    
    
    %convolved function
    kq1 = 0;
    kq2 = 0;
    kq3 = 0;
    
    
    func = @(B1,B2,B3,sigma,t0,k1,k2,k3,x) ((B1.*exp(0.5.*(kq1+k1).*((kq1+k1).*sigma.^2 + 2.*t0 - 2.*x)).*(1 - sqrt(1./sigma.^2).*sigma.*erf(((kq1+k1).*sigma.^2 + t0 - x)./(sqrt(2).*sigma))))./(2.*sqrt(1./sigma.^2).*sigma) + (B2.*exp(0.5.*(kq2+k2).*((kq2+k2).*sigma.^2 + 2.*t0 - 2.*x)).*(1 - sqrt(1./sigma.^2).*sigma.*erf(((kq2+k2).*sigma.^2 + t0 - x)./(sqrt(2).*sigma))))./(2.*sqrt(1./sigma.^2).*sigma) + (B3.*exp(0.5.*(kq3+k3).*((kq3+k3).*sigma.^2 + 2.*t0 - 2.*x)).*(1 - sqrt(1./sigma.^2).*sigma.*erf(((kq3+k3).*sigma.^2 + t0 - x)./(sqrt(2).*sigma))))./(2.*sqrt(1./sigma.^2).*sigma));
    
    %curve fit
    [curvefit,gof] = fit(x,y,func, 'upper', upperbounds, 'lower', lowerbounds, 'StartPoint', startguess,  'Robust', 'LAR','MaxFunEvals',500000,'TolX',10E-8)

    %retrieve fit coefficients
    coeffs = coeffvalues(curvefit);
    B1 = coeffs(1);
    B2 = coeffs(2);
    B3 = coeffs(3);
    sigma = coeffs(4);
    t0 = coeffs(5);
    k1 = coeffs(6);
    k2 = coeffs(7);
    k3 = coeffs(8);

    p = confint(curvefit);
    p = p(2,:);
    n = coeffs-p;
    n(6) = 1./(coeffs(6))-1./p(6);
    n(7) = 1./(coeffs(7))-1./p(7);
    n(8) = 1./(coeffs(8))-1./p(8);

    coeffs(6) = 1./coeffs(6);
    coeffs(7) = 1./coeffs(7);
    coeffs(8) = 1./coeffs(8);
    
     %percentage amplitude
     b1 = (B1)./(B1+B2+B3)
     b2 = (B2)./(B1+B2+B3)
     b3 = (B3)./(B1+B2+B3)
   
    %values
    vals = {'B1','B2','B3','sigma','t0','tau1','tau2','tau3'};

%results table
results = [coeffs' n']

subplot(2,1,1)
hold on
plot(x,func(B1,B2,B3,sigma,t0,k1,k2,k3,x),'b','DisplayName','Design 4');axis([-1 15 0 1.1])
subplot(2,1,2)
hold on
plot(x,func(B1,B2,B3,sigma,t0,k1,k2,k3,x),'b');axis([-1 15 1e-5 1.1])
    end
end
% 
% T = table(x,y);
% writetable(T,'C:\Users\an19697\OneDrive - University of Bristol\Desktop\Flavins\TCSPC\D4f_52k.txt');
