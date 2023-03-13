%%  Kepler optimization algorithm (KOA) source codes demo 1.0               %
%                                                                         %
%  Developed in MATLAB R2019a                                    %
%                                                                         %
%  Author and programmer: Reda Mohamed (E-mail: redamoh@zu.edu.eg) & Mohamed Abdel-Basset (E-mail: mohamedbasset@ieee.org)                              %
%                                                                         %
%   Main paper: Abdel-Basset, M., Mohamed, R.                                    %
%               Kepler optimization algorithm: A new metaheuristic algorithm inspired by Kepler’s laws of planetary motion                         %
%               Knowledge-Based Systems
%               DOI: doi.org/10.1016/j.knosys.2023.110454                 %
%_________________________________________________________________________%
%%
% This function initialize the first population of search agents
function Positions=initialization(SearchAgents_no,dim,ub,lb)

Boundary_no= length(ub); % numnber of boundaries

% If the boundaries of all variables are equal and user enter a signle
% number for both ub and lb
if Boundary_no==1
     Positions=rand(SearchAgents_no,dim).*(ub-lb)+lb;
end

% If each variable has a different lb and ub
if Boundary_no>1
    for i=1:dim
        ub_i=ub(i);
        lb_i=lb(i);
         Positions(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;      
    end
end