clear; close all;

% This script is used to fit a linear incompressible material model to test
% data. Material test data for two different mixing ratios of Sylgard
% CY52-276 are available. The test data for mixing ratios A:B of 9:10 and
% 10:9 can be selected below in the section "Load data". This script
% converts the nominal stress given in the material tests to cauchy stress
% and does a combined least square fit of both the uniaxial as well as the
% equibiaxial test data to a linear model. The parameter "weighting" in the
% section "Coupled least square fit" can be adjust to give relatively less
% weight to the equibiaxial data. The important extracted material
% parameter is the Young's modulus E given in MPa.


%% Load data

load('StrainStressTest_9_to_10.mat')

nu = 0.49; % Poisson's Ratio

%% Calculate cauchy stress from input which is nominal stress
u.strain = uniaxial(:,1);
u.lambda1 = u.strain + 1;
u.lambda2 = 1 ./ sqrt(u.lambda1);
u.lambda3 = u.lambda2;

e.strain = equibiaxial(:,1);
e.lambda1 = equibiaxial(:,1) + 1;
e.lambda2 = e.lambda1;
e.lambda3 = 1./(e.lambda1.*e.lambda2);

for i = 1:length(u.lambda1)
    u.nominal_stress_tensor = zeros(3);
    u.nominal_stress_tensor(1,1) = uniaxial(i,2);
    u.nominal_stress(i,1) = u.nominal_stress_tensor(1,1);
        
    u.F(:,:,i) = diag([u.lambda1(i),u.lambda2(i),u.lambda3(i)]);
    u.cauchy_stress_tensor = u.F(:,:,i)*u.nominal_stress_tensor;
    u.cauchy_stress(i,1) = u.cauchy_stress_tensor(1,1);
    
end
for i = 1:length(e.lambda1)
    e.nominal_stress_tensor = zeros(3);
    e.nominal_stress_tensor(1,1) = equibiaxial(i,2);
    e.nominal_stress_tensor(2,2) = equibiaxial(i,2);
    
    e.nominal_stress(i,1) = e.nominal_stress_tensor(1,1);

    e.F(:,:,i) = diag([e.lambda1(i),e.lambda2(i),e.lambda3(i)]);
    e.cauchy_stress_tensor = e.F(:,:,i)*e.nominal_stress_tensor;
    e.cauchy_stress(i,1) = e.cauchy_stress_tensor(1,1);
end

fx_u = u.cauchy_stress;
x_u = u.strain;
fx_e = e.cauchy_stress;
x_e = e.strain;

%% Coupled least square fit for uniaxial and equibiaxial
weighting = 1; % [0,1] This factor give more weight to uniaxial if lowered

fun1 = @(a)(sum((fx_u-a*x_u)).^2 +...
    sum((fx_e(1:round(numel(fx_e)*weighting))-a*x_e(1:round(numel(fx_e)*weighting))/(1 - nu)).^2));
options = optimset('Display','iter','TolFun',1e-12,'TolX',1e-12);
E = fminsearch(fun1,1,options); % E-Modulus [MPa]

%% Plotting results Cauchy Stress
figure
subplot(1,2,1)
plot(u.strain,u.cauchy_stress,'LineWidth',3,'Color',[41, 128, 185]./255)
ylabel("Cauchy Stress [MPa]")
xlabel("Strain [-]")
title("Uniaxial")
hold on
plot(u.strain,E*u.strain,'LineWidth',1.5,'Color',[39, 174, 96]./255)
legend({'Data','Fit'},'Location','northwest')
ylim([0,max(e.cauchy_stress)])

subplot(1,2,2)
plot(e.strain,e.cauchy_stress,'LineWidth',3,'Color',[41, 128, 185]./255)
ylabel("Cauchy Stress [MPa]")
xlabel("Strain [-]")
title("Equibiaxial")
hold on
plot(e.strain,E*e.strain/(1-nu),'LineWidth',1.5,'Color',[39, 174, 96]./255)

%% Convert cauchy to nominal
u.fit_cauchy_stress = E*u.strain;
e.fit_cauchy_stress = E*e.strain/(1-nu);


for i = 1:numel(u.fit_cauchy_stress)
    f1 = u.F(:,:,i)^-1;
    u.fit_nominal_stress(i) = u.fit_cauchy_stress(i) * f1(1,1);
    
    f1 = e.F(:,:,i)^-1;
    e.fit_nominal_stress(i) = e.fit_cauchy_stress(i) * f1(1,1);
end

%% Plotting results Nominal Stress
figure
subplot(1,2,1)
plot(u.lambda1,u.nominal_stress,'LineWidth',3,'Color',[41, 128, 185]./255)
hold on
plot(u.lambda1,u.fit_nominal_stress,'LineWidth',1.5,'Color',[39, 174, 96]./255)
legend({'Data','Fit'},'Location','northwest')
ylabel("Nominal Stress [MPa]")
xlabel("Stretch [-]")
title("Uniaxial")
ylim([0,max(e.nominal_stress)])

subplot(1,2,2)
plot(e.lambda1,e.nominal_stress,'LineWidth',3,'Color',[41, 128, 185]./255)
hold on
plot(e.lambda1,e.fit_nominal_stress,'LineWidth',1.5,'Color',[39, 174, 96]./255)
ylabel("Nominal Stress [MPa]")
xlabel("Stretch [-]")
title("Equibiaxial")