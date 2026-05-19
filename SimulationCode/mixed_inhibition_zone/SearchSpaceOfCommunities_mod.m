% This code explores the space of all communities for a given N, M and saves the communities which have stable internal fixed points. The space of all communities can be % efficiently traversed by the algorithm outlined in Appendix \ref{appendixA} of our manuscript. The results are saved as a cell array 'communities_with_stable_coexistence' which contains a community index k once for every fixed point this algorithm finds. unique(communities_with_stable_coexistence) will give a list of communities that have at least one stable internal fixed point.

addpath("/Users/prateekverma/Documents/GitHub/antibiotics_biodiversity/SimulationCode/mixed_inhibition_zone/functions")

N = 3; % number of species
M = 2; % number of antibiotics
g = repelem(1,N); % base fitness of each species

% metabolic costs of: 
c_res = repelem(1/20,M);  %intrinsic resistance
c_deg = 2.1*c_res ; %resistance via degradation
c_prod = 3*c_res;  %production

K_D = 10;      % degrading area
K_P = 15;      % killing area

simplex_sampling_size = 200; % number of different initial conditions to supply to solver
timesteps = 600; % number of timesteps to evolve the system given some initial condition. used later. 

% Given a sequence of N*M phenotypes of the N species with respect to the M
% antibiotics, generate the matrices P, S, D, R for each community. 

% % Example: Below are the matrices corresponding to the cyclical interactions in Kelsic et al
% P = [1 0 0; 0 1 0; 0 0 1];
% S = [0 0 1; 1 0 0; 0 1 0];
% D = [0 1 0; 0 0 1; 1 0 0];
% R = [0 0 0; 0 0 0; 0 0 0];

phenotype_sequences = NonRedundantCommunities(N,M); % find set of nonredundant communities

% for a few cases, the folder NRcommunities already contains the NR set
% of communities. It can be called as follows.
% load("NRcommunities/phenseq3,2.mat","phenotype_sequences")
% if line 30 is uncommented, line 26 is not required

% phenotype_sequences = NonRedundantCommunities(N,M);
filename = "NRcommunities/phenseq" + N + "," + M + ".mat";
load(filename,"phenotype_sequences")
cell_counter = 1;
communities_with_stable_coexistence = {};

% Preallocate results
results = struct('community_id', {}, 'stable', {}, 'PSD_fraction', {}, ...
                 'cyclicity_degree', {}, 'motif_balance_CV', {},'robustness', {});


for community_num = 1:size(phenotype_sequences,2) %iterating over communities i.e, different phenotype matrices P, S, D, R

    % Use the phenotype sequence from the NR set to generate the matrices P,S,D,R
    
    sequence = phenotype_sequences(:,community_num); % string of length N*M with entries in {1,2,3,4}
    P = zeros(N,M);
    S = zeros(N,M);
    D = zeros(N,M);
    R = zeros(N,M);
    
    for matrix_index = 1:N*M % populating the matrices P, S, D, R
        % matrix_index will determine the matrix index (i,j) 
        % where we record the phenotype of strain i with respect to antibiotic j
         
        phenotype = sequence(matrix_index,1);
        colnum = mod(matrix_index,M) + (M - mod(matrix_index,M))*isint(matrix_index/M);
        rownum = 1 + (matrix_index - modifiedmod(matrix_index,M))/M;
        if phenotype == 1
            P(rownum, colnum) = 1;
        elseif phenotype == 2
            S(rownum, colnum) = 1;
        elseif phenotype ==3
            D(rownum, colnum) = 1;
        else 
            R(rownum, colnum) = 1;
        end    
    end

    % === Check PSD motif coverage for each antibiotic ===
    psd_flags = false(1, M);
    for j = 1:M
        psd_flags(j) = any(P(:,j)) && any(S(:,j)) && any(D(:,j));
    end
    PSD_fraction = sum(psd_flags) / M;
    
    if PSD_fraction < 1
        fprintf('Community %d: PSD coverage = %.2f. Skipped.\n', community_num, PSD_fraction);
        continue
    end
    
    % Input:
    % P, S, D: N x M binary matrices (N strains, M antibiotics)
    % Each PSD motif is defined by (producer, sensitive, degrader) for a given antibiotic
    
    % Step 1: Extract all PSD motifs
    motif_list = []; % Each row: [strain_P, strain_S, strain_D, antibiotic_index]

    for j = 1:M
        for i = 1:N
            if P(i,j)
                for k = 1:N
                    if S(k,j)
                        for l = 1:N
                            if D(l,j)
                                motif_list = [motif_list; i, k, l, j];
                            end
                        end
                    end
                end
            end
        end
    end
    
    % Step 2: Build motif interaction graph
    % Each motif is a node; link if they share a strain and involve different antibiotics
    num_motifs = size(motif_list,1);
    adj_matrix = zeros(num_motifs); % motif-to-motif adjacency
    
    for m1 = 1:num_motifs
        for m2 = m1+1:num_motifs
            strain_set_1 = motif_list(m1, 1:3);
            strain_set_2 = motif_list(m2, 1:3);
            antibiotic_1 = motif_list(m1, 4);
            antibiotic_2 = motif_list(m2, 4);
    
            shared_strain = any(ismember(strain_set_1, strain_set_2));
            different_antibiotic = antibiotic_1 ~= antibiotic_2;
    
            if shared_strain && different_antibiotic
                adj_matrix(m1, m2) = 1;
                adj_matrix(m2, m1) = 1; % undirected link
            end
        end
    end

    
    % Step 3: Find connected components in motif graph
    G = graph(adj_matrix);
    component_ids = conncomp(G); % returns vector of component IDs for each motif
    
    % Step 4: For each component, record which antibiotics are represented
    linked_antibiotics = [];
    for c = unique(component_ids)
        motifs_in_component = motif_list(component_ids == c, :);
        antibiotics_in_component = unique(motifs_in_component(:,4));
        linked_antibiotics = [linked_antibiotics; antibiotics_in_component];
    end
    
    linked_antibiotics = unique(linked_antibiotics); % avoid duplication
    
    % Step 5: Compute C'
    C_prime = length(linked_antibiotics) / M;

    cyclicity_degree = C_prime; % Cyclicity metric

    motif_balance_CV = NaN;  % Symmetry metric

    % First we find the fixed points of the dynamical system for the mixed-inhibition zone model.
    fixedpts = zeros(N,1);
    i=1; %to increase the size of the array fixedpts
    
    X = sym('X', [1 N], 'real');    
    modeleqn = symmodel(X,N,M,P,S,D,R,K_D,K_P,g,c_prod,c_deg,c_res)-X;
    eqns = modeleqn == repelem(0,N);   % equations to solve
    
    for k = 1:simplex_sampling_size    % starting solver at different initial conditions 
        solution = vpasolve(eqns, X,'Random', true);   % note: vpa might give multiple solutions
        fixedptarray = round(struct2array(solution) .* 1000)/1000; % each column is a fixed point
        for pt = fixedptarray
            if is_internalfixedpt(pt) == 1 % record only "internal" fixed points i.e. those in the interior of the domain
                fixedpts(:,i) = pt;
                i=i+1;
            end  
        end
    end  
       
    % Now we linearize about each fixed point and analyze stability. 
    
    if isequal(fixedpts,repelem(0,N)')
        fprintf('Community number %d: No internal fixed points. \n', community_num);
        continue
    else
        disp(community_num)
        fixedpts = unique(fixedpts', 'rows')' % neglect fixed points that have been repeated
    end
        
    Jacobian = jacobian(symmodel(X,N,M,P,S,D,R,K_D,K_P,g,c_prod,c_deg,c_res), X); 
    
    for ptcounter = 1:size(fixedpts,2) %iterating over fixed points of the system
        pt = fixedpts(:,ptcounter);
        if isequal(pt,repelem(0,N)') 
            continue
        end
        J_atspecificpoint = subs(Jacobian, X, pt');
        J_numeric = double(J_atspecificpoint);
        if log(max(abs(eig(J_numeric))))<0
            phenotype_matrix = zeros(N,M);
            % to store the interaction profile of this community, we construct an
            % N x M matrix with entries in {1,2,3,4} that combines P, S, D, R
            for row_index = 1:N
                phenotype_matrix(row_index,:) = sequence( 1+(row_index-1)*M : row_index*M,1);
            end
            communities_with_stable_coexistence{cell_counter}=community_num; % record community number
            cell_counter = cell_counter + 1;
            eigvals = eig(J_numeric);
            robustness = 1 - max(abs(eigvals));

            % === Visualize the interaction graph for this stable community ===
            % === Custom Interaction Plot with Annotations ===
            theta = linspace(0, 2*pi, N+1); theta(end) = [];
            x = cos(theta); y = sin(theta);
            strain_labels = arrayfun(@(i) ['Strain ' num2str(i)], 1:N, 'UniformOutput', false);

            figure('Name', ['Community ' num2str(community_num) ' - Interaction Graph']);
            hold on
            axis equal; axis off

            % Plot strain nodes
            for i = 1:N
                plot(x(i), y(i), 'ko', 'MarkerSize', 20, 'LineWidth', 2, 'MarkerFaceColor', 'w');
                text(x(i), y(i), strain_labels{i}, 'FontWeight', 'bold', ...
                    'Color', 'k', 'HorizontalAlignment', 'center', 'FontSize', 12);
            end

            % Plot solid P→S arrows and annotate with Ai
            for j = 1:M
                for p = find(P(:,j))'
                    for s = find(S(:,j))'
                        if p ~= s
                            % Coordinates
                            x_p = x(p); y_p = y(p);
                            x_s = x(s); y_s = y(s);
                            % Direction vector
                            dx = x_s - x_p; dy = y_s - y_p;

                            % Shorten arrow slightly to not hit node centers
                            shrink = 0.1;  % tune if needed
                            L = sqrt(dx^2 + dy^2);
                            dx_adj = (1 - shrink) * dx;
                            dy_adj = (1 - shrink) * dy;

                            % Draw solid arrow
                            quiver(x_p, y_p, dx_adj, dy_adj, 0, 'k', ...
                                'LineWidth', 2, 'MaxHeadSize', 0.1);

                            % Label with A1, A2...
                            mx = (x_p + x_s) / 2;
                            my = (y_p + y_s) / 2;
                            text(mx, my, ['A' num2str(j)], ...
                                'FontSize', 10, 'Color', 'k', ...
                                'FontWeight', 'normal', 'HorizontalAlignment', 'center');

                            % Plot D → midpoint(P→S) as dashed
                            for d = find(D(:,j))'
                                if d ~= p && d ~= s
                                    xd = x(d); yd = y(d);
                                    quiver(xd, yd, mx - xd, my - yd, 0, ...
                                        'LineStyle', '--', 'Color', [0.3 0.3 0.3], ...
                                        'LineWidth', 1.5, 'MaxHeadSize', 0.1);
                                end
                            end
                        end
                    end
                end
            end

            title(['Community ', num2str(community_num), ' Interaction Graph']);
            hold off

            results(1+end) = struct( ...
                'community_id', community_num, ...
                'stable', true, ...
                'PSD_fraction', PSD_fraction, ...
                'cyclicity_degree', cyclicity_degree, ...
                'motif_balance_CV', motif_balance_CV, ...
                'robustness', robustness...
            );
        end
    end

end

% Now we save all communities which have stable fixed points. with multiplicity - 
% if a community has multiple stable fixed points, then that community shows up
% once for every fixed point. 

flag = true; % do you want to save?

if flag == true
    savename = "Results/Search-Space/" + strcat('N=',num2str(N),'_M=', num2str(M)) + ".mat"
    save(savename,'communities_with_stable_coexistence')
end

% Convert results to table
T = struct2table(results);
disp(T)

% Optionally save
if flag
    result_table_name = "Results/Table-N" + N + "-M" + M + ".mat";
    save(result_table_name, 'T')
end



