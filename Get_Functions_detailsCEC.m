%  Kepler optimization algorithm (KOA) source codes demo 1.0               %
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
% This function containts full information and implementations of the benchmark 
% functions in Table 1, Table 2, and Table 3 in the paper

% lb is the lower bound: lb=[lb_1,lb_2,...,lb_d]
% up is the uppper bound: ub=[ub_1,ub_2,...,ub_d]
% dim is the number of variables (dimension of the problem)

function [lb,ub,dim] = Get_Functions_detailsCEC(F)


switch F
    case 1
        dim=10;
        lb=-100*ones(1,dim);
        ub=100*ones(1,dim);
        
    case 2
       dim=10;
        lb=-100*ones(1,dim);
        ub=100*ones(1,dim);
        
    case 3
      dim=10;
        lb=-100*ones(1,dim);
        ub=100*ones(1,dim);
        
    case 4
        dim=10;
        lb=-100*ones(1,dim);;
        ub=100*ones(1,dim);;
        
    case 5
         dim=10;

        lb=-100*ones(1,dim);;
        ub=100*ones(1,dim);;

        
    case 6
         dim=10;

        lb=-100*ones(1,dim);;
        ub=100*ones(1,dim);;

    case 7
        dim=10;

        lb=-100*ones(1,dim);;
        ub=100*ones(1,dim);;
 
    case 8
       dim=10;

        lb=-100*ones(1,dim);;
        ub=100*ones(1,dim);;
  
    case 9
 dim=10;

        lb=-100*ones(1,dim);;
        ub=100*ones(1,dim);;
        
    case 10
         dim=10;

        lb=-100*ones(1,dim);;
        ub=100*ones(1,dim);;
        
    case 11
       dim=10;
        lb=-600*ones(1,dim);;
        ub=600*ones(1,dim);;
        
    case 12
        dim=10;

        lb=-100*ones(1,dim);;
        ub=100*ones(1,dim);;
        
    case 13
        dim=10;

        lb=-100*ones(1,dim);;
        ub=100*ones(1,dim);;
        
    case 14
         dim=10;

        lb=-100*ones(1,dim);;
        ub=100*ones(1,dim);;
        
    case 15
         dim=10;

        lb=-100*ones(1,dim);;
        ub=100*ones(1,dim);;
        
    case 16
         dim=10;

        lb=-100*ones(1,dim);;
        ub=100*ones(1,dim);;
        
    case 17
         dim=10;

        lb=-100*ones(1,dim);;
        ub=100*ones(1,dim);;
        
    case 18
         dim=10;

        lb=-100*ones(1,dim);;
        ub=100*ones(1,dim);;
        
    case 19
         dim=10;

        lb=-100*ones(1,dim);;
        ub=100*ones(1,dim);;
        
        
    case 20
         dim=10;

        lb=-100*ones(1,dim);;
        ub=100*ones(1,dim);;
             
        
    case 21
         dim=10;

        lb=-100*ones(1,dim);;
        ub=100*ones(1,dim);;
            
        
    case 22
         dim=10;

        lb=-100*ones(1,dim);;
        ub=100*ones(1,dim);;
            
        
    case 23
         dim=10;

        lb=-100*ones(1,dim);;
        ub=100*ones(1,dim);;
        
        
  case 24
      dim=10;
        lb=-100*ones(1,dim);;
        ub=100*ones(1,dim);;
        
        
   case 25
         dim=10;

        lb=-100*ones(1,dim);;
        ub=100*ones(1,dim);;
   case 26
        dim=10;

        lb=-100*ones(1,dim);;
        ub=100*ones(1,dim);;
   case 27
        dim=10;

        lb=-100*ones(1,dim);;
        ub=100*ones(1,dim);;
   case 28
         dim=10;

        lb=-100*ones(1,dim);;
        ub=100*ones(1,dim);;
   case 29
         dim=10;

        lb=-100*ones(1,dim);;
        ub=100*ones(1,dim);;
   case 30
         dim=10;

        lb=-100*ones(1,dim);;
        ub=100*ones(1,dim);;
 
 end
end
