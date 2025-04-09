%% Main Community Detection x noisy signal (real or synthetic)
clear, close all
clc
%rng('default')
%s = rng;
addpath(genpath('.\CMRF_EI'))
addpath(genpath( '.\mycommtable_filefolders'))% information theory and functions for CD with classical methods
PC_tiziana = 1; % Camilla da cambiare in 0 se è il tuo PC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parameter settings

% Definition of parameters in params for data generation:
params.generate_trials_sameGamma = 1; %Struttura comunità grafi (vettore)
params.generate_trials_sameA = 0; %Posso cambiarlo con 1
params.isnoisy = 0;
params.isreal = 0; %0

% Graph parameters
% Definition of graph parameters and in G_template for data and graph structure generation
% ---> Attention: if params.isreal = 1, these parameters in G_template are
% useless
G_template.N = 50; %Nodi 100 per sintetici %59
G_template.P = 4; %Comunità
G_template.M = 1; %Layer
G_template.Nruns = 100; %Numero di run %2 % 50 100
G_template.visualize=0; %Visualizzazione %0
G_template.save=0; %salvataggio %0
G_template.snr = 0; %Rapporto segnale-rumore %50 %100 %15

n_iter=20;
G_template.prob=0.95;  %0.99 %probabilità per generazione grafo: prob per Erdos-Renyi, intraprob=1-interprob per modular
G_template.layout='modular';  %tipo di grafo, 'modular' 'erdos_renyi'  'minnesota'
G_template.var_expr='random';    %tipo di segnale da generare, 'random' 'smooth' 'gmrf'
%G_template.var_signoise=0; %varianza rumore
G_template.var_param=0.;

% Graph Learning parameter setting
G_template.selected_alpha=(-0.99.^(1:n_iter)+1)*0.75+0.25;%peso exp square distance w.r.t. comm distance 0.25*ones(n_iter,1);%
G_template.theta_fissato=0.5;
if G_template.prob <=0.5
    G_template.potenziale='repulsivo';
else
    G_template.potenziale='attrattivo';
end
G_template.Markov_pesi_binari=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data generation
params.PC_tiziana = PC_tiziana;
if params.isreal == 1
    params.subj = 'g';
    params.fband = 'beta';
end

params.class='rest';
[array_of_G_class1] = datageneration_05(G_template,params);

% Plot of the original graph

figure ('Name', strcat('Original Graph'))
itrail = 2;%14
G_run =array_of_G_class1{1,itrail};
G_plot = graph(array_of_G_class1{1,itrail}.A_lay);
%h = plot(G_plot,'XData',G_run.coords(:,1),'YData',G_run.coords(:,2), 'EdgeColor',[0.3, 0.3 0.3],'MarkerSize', 10, 'NodeLabel', {});
%h = plot(G_plot,'XData',G_run.coords(:,1),'YData',G_run.coords(:,2), 'EdgeColor',[0.3, 0.3 0.3],'MarkerSize', 10,'NodeLabel', G_run.lab);
h = plot(G_plot,'XData',G_run.coords(:,1),'YData',G_run.coords(:,2), 'EdgeColor',[0.3, 0.3 0.3],'MarkerSize', 10);
colormap('turbo')
box off
axis off
%h.NodeCData = G_run.savedCommunity; %rimetti per




%% Eigenvector analysis
% Define the S_matrix with 6 rows and 5 columns
 % S_matrix = [1500, 1400, 1300, 1200, 1100;
 %     1000, 900, 800, 700, 600;
 %     800, 700, 600, 500, 400;
 %     650, 590, 530, 470, 350;
 %     500, 450, 400, 350, 300;
 %     300, 250, 200, 150, 100];
 % S_matrix = [1500, 1000, 800;
  %   1000, 700, 500;
  %   800, 500, 300;
  %   650, 350, 250;
  %   500, 300, 100;
  %   400, 150, 75]; 
  % 
  % S_matrix = [1500, 1000, 700;
  %    1000, 700, 500;
  %    800, 400, 200;
  %    650, 350, 150;
  %    500, 200, 100;
  %    300, 100, 75]; %ok
  
  S_matrix = [%1500, 800, 300;
     5000, 5000, 5000, 5000, 1000, 500, 100; 
     5000, 5000, 5000, 5000, 1000, 500, 50;
     2000, 2000, 2000, 2000, 1000, 500, 30;
     2000, 2000, 2000, 2000, 1000, 500, 10]
    % 500, 10, 500;
     %325, 5, 200]
     %250, 150, 50;
     %160, 100, 30;
     %80, 50, 15];
     %300, 100, 75]; %ok

   % S_matrix = [
   %   800, 400, 200;
   %   400, 200, 100;
   %   200, 80, 30;
   %   150, 60, 10
   %   ]; % ok per rete piccola

  % S_matrix = [1500,  1100, 900;
  %    1000,  500, 100;
  %    800,  300,  80;
  %    650,  200,  90;
  %    500,  100,  80;
  %    300,  100,  75];
    % S_matrix = [1500,  1100, 900;
    %  1000,  800, 500;
    %  800,  300,  100;
    %  650,  200,  90;
    %  300,  100,  80;
    %  100,  75,  40];
 % S_matrix = [1500, 1500, 1500, 1500, 1500;
 %     1000, 1000, 1000, 1000, 1000;
 %     750, 750, 750, 750, 750;
 %     250, 250, 250, 250, 250;
 %     100, 100, 100, 100, 100;
 %     75, 75, 75, 75, 75];

% Loop over the rows of S_matrix (6 rows)
for s = 1:size(S_matrix, 1)
    for g = 1:G_template.Nruns
        G_run = array_of_G_class1{1, g};
        A = G_run.A_lay;
        non_zero_mask = A ~= 0;

        % Loop over each layer and apply the transformation and perturbation in one step
        for k = 1:size(S_matrix, 2)
            S_value = S_matrix(s, k); % Extract the specific S value for transformation

            % Define the sigmoid-like transformation centered at 0.5
            transformation = @(x) 1 ./ (1 + exp(-S_value * (x - 0.5)));

            % Apply the transformation only to non-zero elements for the current layer
            A_transformed = zeros(size(A));
            A_transformed(non_zero_mask) = transformation(A(non_zero_mask));
            %A = A_transformed;

            theta = 0.5;

            % Binarize the adjacency matrix
            B = A_transformed >= theta;
            stability_edge = abs(A_transformed - theta);
            B_L = diag(sum(B, 1)) - B;

            [U, Lambda] = eig(double(B_L)); % Eigenvectors of nominal binary graph

            % Compute probability P_m
            P_m = zeros(size(A_transformed));

            matrix_overthreshold_links = (A_transformed>=theta);
            matrix_underthreshold_links = (A_transformed<theta);

            P_m(matrix_overthreshold_links) = 1-A_transformed(matrix_overthreshold_links); % forse mettere A-1?
            P_m(matrix_underthreshold_links) = A_transformed(matrix_underthreshold_links);
            
            gamma_vect_k(:) = computegamma(U,P_m); 
            
            Niter = 100;

            % Loop over iterations for perturbations
            for iiter = 1:Niter
                % Preallocate arrays for storing results
                random_var = 0.5 * rand(G_run.N);
                random_var = triu(random_var) + triu(random_var, 1)';

                community_matrix = G_run.savedCommunity(:);
                inter_community_mask = community_matrix ~= community_matrix';
                
                random_var = random_var .* inter_community_mask;

                % Determine which edges to perturb based on the stability edge
                modified_edge = random_var > stability_edge;

                num_modifiededge(s,k,iiter,g) = sum(modified_edge(:));


                % Perturb the Laplacian matrix
                B_pert = B_L;
                indices = modified_edge == 1;
                B_pert(indices) = 1 - B_L(indices);

                % Calculate the perturbed eigenvalues and eigenvectors
                [U_pert(:,:,iiter), Lambda_pert(:,:,iiter)] = eig(double(B_pert));

                

                binary_L_pert(:,:,iiter) = B_pert;
                
                gamma_vect_all(:,iiter,k) = gamma_vect_k;
             end

            U_pert_all(:,:,:,k) = U_pert;
            Lambda_pert_all(:,:,:,k) = Lambda_pert;
            binary_L_pert_all(:,:,:,k) = binary_L_pert;

            U_nominal_allk(:,:,k) = U;
            Lambda_nominal_allk(:,:,k) = Lambda;

        end

        % 

        mean_binary_L_pert_alliter = mean(binary_L_pert_all,4); % ci sono tutte le iiter
        U_meanU_pert = mean(U_pert_all,4);
        Lambda_meanU_pert = mean(Lambda_pert_all,4);

        choice_nom_graph = 1;
        for iiter = 1:Niter
            
            alpha_2_k = gamma_vect_all(2,iiter,:).^(-1)/(sum(gamma_vect_all(2,iiter,:).^(-1)));
            alpha_4_k = gamma_vect_all(4,iiter,:).^(-1)/(sum(gamma_vect_all(4,iiter,:).^(-1)));
            Lambda_2_smart = squeeze(alpha_2_k)'*squeeze(Lambda_pert_all(2,2,iiter,:));

            %U_pert_all(N,N,Niter,Klayers)
            U_smart = zeros(G_template.N,G_template.P);
            Lambda_smart = zeros(G_template.N,1);

            K = size(S_matrix, 2);
            P = G_template.P;
            N = G_template.N;
            
            for myp = 1:P
                alpha_myk = gamma_vect_all(myp,iiter,:).^(-1)/(sum(gamma_vect_all(myp,iiter,:).^(-1)));
                for myk = 1:K
                    U_smart(:,myp) = U_smart(:,myp) + U_pert_all(:,myp,iiter,myk)*alpha_myk(myk);
                end
                U_smart(:,myp) = U_smart(:,myp)/sqrt(sum(U_smart(:,myp).^2));
                U_Unorm(:,myp) = U_meanU_pert(:,myp,iiter)/sqrt(sum(U_meanU_pert(:,myp,iiter).^2));
                
                for myn = 1:N
                    %Lambda_smart(myn) = squeeze(alpha_myk)'*squeeze(Lambda_pert_all(myk,myn,iiter,:)); % mi sembra da sistemare questo
                    Lambda_smart(myn) = squeeze(alpha_myk)'*squeeze(Lambda_pert_all(myn,myn,iiter,:));
                end
            

            end


            
            [U_meanL_pert(:,:,iiter), Lambda_meanL_pert(:,:,iiter)] = eig( mean_binary_L_pert_alliter(:, :, iiter));

            %dist_U_nom(iiter) = norm(U_meanU_pert(:,1:G_template.P ,iiter) - U_nominal_allk(:,1:G_template.P,choice_nom_graph), 'fro');
            
            dist_U_nom(iiter) = grasdist (U_Unorm(:,3:end),U_nominal_allk(:,3:P,choice_nom_graph)); 

            %dist_L_nom(iiter) = norm(U_meanL_pert(:,1:G_template.P ,iiter) - U_nominal_allk(:,1:G_template.P,choice_nom_graph), 'fro');
            dist_L_nom(iiter) = grasdist(U_meanL_pert(:,3:P,iiter), U_nominal_allk(:,3:P,choice_nom_graph));

            %dist_Usmart_nom(iiter) = norm(U_smart- U_nominal_allk(:,1:G_template.P,choice_nom_graph), 'fro');
            dist_Usmart_nom(iiter) =grasdist(U_smart(:,3:end), U_nominal_allk(:,3:P,choice_nom_graph));

            NPM_meanUvsnom(:,iiter) = abs(diag(Lambda_meanU_pert(:,:,iiter))-diag(Lambda_nominal_allk(:,:,choice_nom_graph)))./diag(Lambda_nominal_allk(:,:,choice_nom_graph));
            NPM_meanLvsnom(:,iiter) = abs(diag(Lambda_meanL_pert(:,:,iiter))-diag(Lambda_nominal_allk(:,:,choice_nom_graph)))./diag(Lambda_nominal_allk(:,:,choice_nom_graph));
            NPM_meanUsmartvsnom(:,iiter) = abs(Lambda_smart-diag(Lambda_nominal_allk(:,:,choice_nom_graph)))./diag(Lambda_nominal_allk(:,:,choice_nom_graph));


        end

        fro_d_UmeanU_nominal_mean(g,s) = mean(dist_U_nom);
        fro_d_Usmart_nominal_mean(g,s) = mean(dist_Usmart_nom);
        fro_d_UmeanL_nominal_mean(g,s) = mean(dist_L_nom);

        NPM_Lmean_nominal_mean(:,g,s) = mean(NPM_meanLvsnom,2);
        NPM_Umean_nominal_mean(:,g,s) = mean(NPM_meanUvsnom,2);
        NPM_Usmart_nominal_mean(:,g,s) = mean(NPM_meanUsmartvsnom,2);

    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Performance evaluation (different types of average)

num_modifiededge_mean = mean(mean(mean(num_modifiededge,4),3),2);
% eigenvectors
fro_Lmean = mean(fro_d_UmeanL_nominal_mean,1);
fro_Umean = mean(fro_d_UmeanU_nominal_mean,1);
fro_Usmart = mean(fro_d_Usmart_nominal_mean,1);


figure
plot(num_modifiededge_mean,fro_Lmean, '-o','LineWidth',3,'Color', [0.4660 0.6740 0.1880]);
hold on
plot(num_modifiededge_mean, fro_Umean,'-o','LineWidth',3,'Color',  [0.9290 0.6940 0.1250]);
hold on
plot(num_modifiededge_mean, fro_Usmart,'-o','LineWidth',3,'Color', [0.3010 0.7450 0.9330]);

grid on
legend('Mean on L', 'Mean on U','Weighted mean on U');
xlabel('number of perturbed edges','FontName', 'Times New Roman')
ylabel('subspace distance','FontName', 'Times New Roman')
set(gca,'fontsize', 14,'FontName', 'Times New Roman');
hold off

% eigenvalues
NPM_Lmean = squeeze(mean(NPM_Lmean_nominal_mean,2));
NPM_Umean = squeeze(mean(NPM_Umean_nominal_mean,2));
NPM_Usmart = squeeze(mean(NPM_Usmart_nominal_mean,2));

figure
plot(num_modifiededge_mean,sum(NPM_Lmean(3:P,:),1), '--o','LineWidth',3, 'Color', [0.4660 0.6740 0.1880]);
hold on
plot(num_modifiededge_mean, sum(NPM_Umean(3:P,:),1),'--o','LineWidth',3, 'Color',  [0.9290 0.6940 0.1250]);
hold on
plot(num_modifiededge_mean, sum(NPM_Usmart(3:P,:),1),'--o','LineWidth',3, 'Color', [0.3010 0.7450 0.9330]);
hold on

 %plot(num_modifiededge_mean,NPM_Lmean(3,:), '-o','LineWidth',3, 'Color', [0.4660 0.6740 0.1880]);
 %hold on
 %plot(num_modifiededge_mean, NPM_Umean(3,:),'-o','LineWidth',3, 'Color', [0.9290 0.6940 0.1250]);
 %hold on 
 %plot(num_modifiededge_mean, NPM_Usmart(3,:),'-o','LineWidth',3, 'Color', [0.3010 0.7450 0.9330]);
 %hold on

grid on
legend('Mean on L', 'Mean on \lambda ', 'Weighted mean on \lambda');
%legend('Mean on L NPM_2', 'Mean on U NPM_2', 'Smart U NPM_2','Mean on L NPM_4', 'Mean on U NPM_4','Smart U NPM_4');
xlabel('number of perturbed edges','FontName', 'Times New Roman')
ylabel('\Sigma_{i=1}^P NEM_i','FontName', 'Times New Roman')
set(gca,'fontsize', 14,'FontName', 'Times New Roman');
hold off




%% functions
function my_Csum_m = computegamma(U,P_m) %scandisco i layer

N = size(U,1);
for i = 1:N
    diffmat_i = (U(:,i)-U(:,i)');
    diffsqmat_i = diffmat_i.^4.*P_m.*(1-P_m);
    my_Cmat = [diffsqmat_i];
    my_Csum_m(i) = sum(my_Cmat(:));

end
end

function grassmann_distance = grasdist(V,U)
P1 = V * (V' * V) ^ -1 * V';
P2 = U * (U' * U) ^ -1 * U';

% Calcola la matrice risultante della proiezione del primo sottospazio sul secondo
M = P1 * P2;

S = svd(V'*U);
%S = min(max(S, -1), 1);
angles = acos(S);

% Somma dei moduli

% La distanza di Grassmann è la norma (norma 2) degli angoli principali
grassmann_distance = sqrt(sum(angles.^2));

P = size(V,2);
for i=1:P
    dist(i)= sqrt(sum((V(:,i) - U(:,i)).^2));
end
%grassmann_distance = sum(dist);

% theta = subspace(V, U);
%grassmann_distance = sin(theta);
end

