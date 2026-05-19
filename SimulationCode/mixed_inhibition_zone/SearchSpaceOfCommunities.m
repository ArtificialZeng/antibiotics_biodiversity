% This code explores the space of all communities for a given N, M and saves the communities which have stable internal fixed points. The space of all communities can be % efficiently traversed by the algorithm outlined in Appendix \ref{appendixA} of our manuscript. The results are saved as a cell array 'communities_with_stable_coexistence' which contains a community index k once for every fixed point this algorithm finds. unique(communities_with_stable_coexistence) will give a list of communities that have at least one stable internal fixed point.

addpath("../functions/")

N = 3; % number of species
M = 2; % number of antibiotics
g = repelem(1,N); % base fitness of each species

% metabolic costs of: 
c_res = repelem(1/20,M);  %intrinsic resistance
c_deg = 2.1*c_res ; %resistance via degradation
c_prod = 3*c_res;  %production

K_D = 10;      %degrading area
K_P = 15;     %killing area

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

for community_num = 1:size(phenotype_sequences,2) %iterating over communities i.e, different phenotype matrices P, S, D, R

    % Use the phenotype sequence from the NR set to generate the matrices P,S,D,R
    
    sequence = phenotype_sequences(:,community_num); % string of length N*M with entries in {1,2,3,4}
    P = zeros(N,M);
    S = zeros(N,M);
    D = zeros(N,M);
    R = zeros(N,M);
    
    for matrix_index = 1:N*M %populating the matrices P, S, D, R
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
    
    % Uncomment the following block to plot trajectory of the system from 5
    % separate initial conditions uniformly distributed in the simplex. 
    
%     for init = UniformSampleSimplex(N,5)
%         init
%         trajectory = TimeEvolution(init, N,M,P,S,D,R,K_D,K_P,g,c_prod,c_deg,c_res, timesteps);
%         for coordinate=1:N
%             plot(1:timesteps,trajectory(coordinate,:),'DisplayName',num2str(coordinate),'LineWidth',3)
%             hold on
%         end
%         legend
%         title("Community number "+community_num)
%         xlabel('Time')
%         ylabel('X')
%         ylim([0,1])
%         drawnow
%         pause(15)
%         hold off
%     end

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
        if log(max(abs(eig(J_atspecificpoint))))<0
            phenotype_matrix = zeros(N,M);
            % to store the interaction profile of this community, we construct an
            % N x M matrix with entries in {1,2,3,4} that combines P, S, D, R
            for row_index = 1:N
                phenotype_matrix(row_index,:) = sequence( 1+(row_index-1)*M : row_index*M,1);
            end
            communities_with_stable_coexistence{cell_counter}=community_num; % record community number
            cell_counter = cell_counter + 1;
            
            % Uncomment the following block to plot trajectory of the three coordinates after
            % starting near pt (optionally after perturbation from pt). 
            
%             initial = pt; % to perturb, add [[-0.01,0.01] zeros(1,numel(pt))] to pt
%             trajectory = TimeEvolution(initial, N,M,P,S,D,R,K_D,K_P,g,c_prod,c_deg,c_res, timesteps);
%             for coordinate=1:N
%                 plot(1:timesteps,trajectory(coordinate,:),'DisplayName',num2str(coordinate),'LineWidth',3)
%                 hold on
%             end
%             legend
%             title("Community number "+community_num)
%             xlabel('Time')
%             ylabel('X')
%             drawnow
%         else
%             initial = pt;
%             trajectory = TimeEvolution(initial, N,M,P,S,D,R,K_D,K_P,g,c_prod,c_deg,c_res, timesteps);
%             for coordinate=1:N
%                 plot(1:timesteps,trajectory(coordinate,:),'DisplayName',num2str(coordinate), 'LineWidth',3)
%                 hold on
%             end
%             legend
%             title("Community number "+community_num)
%             xlabel('Time')
%             ylabel('X')
%             drawnow
        end
%         pause(10)
%         hold off
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

 



