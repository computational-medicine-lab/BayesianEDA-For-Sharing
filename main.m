% tonic model research
% model without chi
clear all; close all; clc;
warning('off','all');

addpath('dependencies');
data = []; 

tau = [1 4 5];
theta_ub = [6 300 10]'; theta_lb = [1.5 50 0.1]'; theta_mu = [2.77325 100 1]'; theta_var = [(0.521739)^2 (10)^2 2^2];

Fss = 100;
Fsu = 4;
Fs = Fsu; Fsy = Fs;
Ts = 1/Fs;

signal_starting_point = 150;
signal_duration = 200;

%% load data

malesub = [11,26,8,10,20,23,3,4,9,13,17,22,24];
femalesub = [12,15,7,18,21,25,1,2,5,6,14,16,19];

sub_list = [malesub,femalesub];

load('ss_8_14_18.mat');
for subject = sub_list
    close all;
    
    %% low pass filter and downsample the data
    data.subject = subject;
    y = ss(subject).sc(1).y;
    noise_sigma = get_sigma(y, 100, 30);
    data.noise_sigma_max = noise_sigma;
    
    y = LowPassFilter(y, Fss, 2,4096);
    
    y = downsample(y,Fss/Fsy); % downsample to 4 Hz
    
    % take only the deconvolution segment of interest
    y = y(signal_starting_point*Fsy:(signal_starting_point+signal_duration)*Fsy-1); y = y(:)';

    noise_dB = 20*log10(var(y)/noise_sigma^2);


%% initialize the structure with the parameters
    data.GCV_WIN = 100;

    data.model = '3D';
    data.chi = 0.50; % assuming (chi * 100) % of the rate of change of amount of sweat in the sweat duct is because of diffusion

    % set model
    data.SC_model = @(theta__, Ts__)SC_model(theta__, Ts__);

    % set constrains and priors
    data = set_constrains(data);

    % set theta intializer
    data.rand_init = @(data_)rand_init(data_);

    data.epsilon = 1e-5; % perturbation

    data.Ts = Ts; % sampling period
    data.Fsu = Fsu;% sampling freq for u
    data.Fsy = Fsy;% sampling freq for y
    data.Fs = Fs;% original freq Fs = Fsy = Fsu in this case, same as before
    data.signal_duration = signal_duration; % the duration of the signal that needs to be deconvolved
    data.noise_sigma = 1e-4;%noise_sigma / 500; % value of measurement noise
    
    data.irls_itermax = 10; %
    data.iteration_for_init = 10;
    % data.lambda_par = 0;
    data.lambda_max = 1e1; % for IRLS
    data.lambda_init = 1e-7; % for IRLS
    data.max_iter = 1000; % max iter of EM
   
    data.u_th = 0.03; % threshold during sparse recovery
    data.u_th_IRLS = 0.03; % threshold during sparse recovery
    data.u_th_GCV = 0.25; % threshold during sparse recovery
    
    data.conv_criteria = 1; % percent

    data.y_obs = y; %observation
    data.x(:,1) = [0;0;data.y_obs(:,1)]; % intial condition
    number_of_initializations = 1; % number of random initialization for the parameters to account for the non-convexity. No multiple initialization is needed as results converges to the same location.
    data.display_flag = 0; % turn off/on display, 1 means on.
    data.estimate_observation_noise_using_EM = true; % logical flag for estimating 'noise_sigma' or use the predefined value.
    
    %% Perform deconvolution
    tic
    for k = 1:number_of_initializations
        data_results(k) = bayesianEDA(data);
    end

    %% select the best results that minimizes the cost
    cost = Inf;
    for k = 1:number_of_initializations

        %if(data_results(k).convergence_flag == 1)% && data_results(k).booundary_stagnation_flag == 0
            if(data_results(k).cost2 < cost)
                data = data_results(k);
                cost = data_results(k).cost2;
            end
        %end
    end

    toc
    %% draw figures with CI
    [hh1, hh2] = plot_results(data);
    


    SNR = 10*log(var(y)/noise_sigma^2);

%% save workspace
    result_dir  = 'Experimemntal_Data_Results';
    mkdir(result_dir)
    % 
    save([result_dir,'\result',num2str(subject)]);
    
    saveas(hh1,[result_dir,'\S',num2str(subject),'_Deconv','.png']);
    saveas(hh2,[result_dir,'\S',num2str(subject),'_conv','.png']);
   
end



%% model
function [Ad, Bd, Cd, Dd, Cd_tonic, Cd_Phasic] = SC_model(theta, Ts)


%% model 1: 3D
% A = [-1/theta(1) 0 0; +1/theta(1) -1/theta(2) 0; +1/(theta(1)) 0 -1/(theta(3))];
% B = [1;0; 0];  
% % B = [1;0; 0];  
% C = [0 theta(4) theta(5)];  D = 0;
% C_tonic =  C .* [0 0 1] ; C_phasic = C .* [0 1 0];

% model 2: 3D
A = [-1/theta(1) 0 0; +1/theta(1) -1/theta(2) 0; +1/(theta(1)) 0 -1/(theta(3))];
B = [1;0; 0];  
% B = [1;0; 0];  
C = [0 theta(4) theta(5)];  D = 0;
C_tonic =  C .* [0 0 1] ; C_phasic = C .* [0 1 0];

% %% model 3: 3D
% A = [-1/theta(1) 0 0; +1/theta(2) -1/theta(2) 0; +1/(theta(2)) 0 -1/(theta(3))];
% B = [1;0; 0];  
% % B = [1;0; 0];  
% C = [0 theta(4) theta(5)];  D = 0;
% C_tonic =  C .* [0 0 1] ; C_phasic = C .* [0 1 0];
%% model 4: 2D
% A = [-1/theta(1) 0; +1/theta(1) -1/theta(2)];
% B = [1;0];  
% % B = [1;0; 0];  
% C = [1 theta(3)];  D = 0;
% C_tonic =  C .* [0 1] ; C_phasic = C .* [1 0];

%% symbolic calculation for discretization
% tic
% syms t 
% expr = expm(A*t); F = int(expr,0,Ts);
% sysd.A = expm(A*Ts);  sysd.B = double(subs(F))*B;
% sysd.C = C; sysd.D = D;
% toc
%% matlab c2d for discretization
sys = ss(A,B,C,D); 
sysd = c2d(sys, Ts, 'impulse');
% sysd = c2d(sys, Ts, 'tustin');
Ad = sysd.A; Bd = sysd.B; Cd = sysd.C; Dd = sysd.D;

Cd_tonic = C_tonic; Cd_Phasic = C_phasic;
end

%% optimization constraints and probabilistic priors
function data = set_constrains(data)
    if(data.model == '3D')
    %% model 1: 3D 
        data.constrains.C = [-1 0 0 0 0;
            0 -1 0 0 0;
            0 0 -1 0 0;
            2 -1 0 0 0;
            0 15 -1 0 0;];

        data.constrains.b = [-0.2;-0.1;-6;0;0];
        data.constrains.Ceq = [0 0 0 1 0; 0 0 0 0 1];
        data.constrains.beq = [1-data.chi; data.chi];

        data.theta_ub = []';
        data.theta_lb = []';


        data.theta_mu = [0.65 2.77325 120 1-data.chi data.chi]';
        data.theta_g = data.theta_mu;
        data.theta_var = [(0.212443)^2 (0.521739)^2 80^2 3^2 3^2]';

        data.lambda_par = [3e1, 1e1, 2e-1, 0, 0]'; %for the theoratical maximization


    elseif(data.model == '2D')
        %% model 1: 2D 
            data.constrains.C = [-1 0 0; 0 -1 0; 0 0 -1];
            data.constrains.b = [0; 0; 0];
            data.constrains.Ceq = [];
            data.constrains.beq = [];

            data.theta_ub = [6 300 10]';
            data.theta_lb = [1.5 50 0.1]';            
            data.theta_mu = [2.77325 100 1];
            data.theta_var = [(0.521739)^2 (10)^2 2^2];
            data.lambda_par = [1e6, 1e6, 0]';
    end
end

%% function for parameter random initialization
function theta = rand_init(data)
    rng('shuffle')
    random_init_flag = 1;
    if(data.model == '3D')
        while(random_init_flag)
            theta(1) =   0.1 + rand(1,1); %[0.70 3 100 0.4];
            theta(2) =  theta(1) * 8; 
            theta(3) = theta(2) * 15; 
            theta(5) = data.chi;
            theta(4) = 1-data.chi;
            theta = theta(:);

            [random_init_flag,~,~] = (check_feasibility(theta,data.constrains.C,data.constrains.b,data.constrains.Ceq,data.constrains.beq));
            random_init_flag = 1-random_init_flag;
            if random_init_flag == 1
                disp('Infisible Initialization')
            end
        end
    elseif(data.model == '2D')
            theta(1) =   0.1 + rand(1,1); %[0.70 3 100 0.4];
            theta(2) = 0.1 + 100*rand(1,1);
            theta(3) = 0.1 + 3*rand(1,1); 
    end


end