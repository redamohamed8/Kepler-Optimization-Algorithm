%_________________________________________________________________________%
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


% The Kepler Optimization Algorithm
function [Sun_Score,Sun_Pos,Convergence_curve]=KOA(SearchAgents_no,Tmax,ub,lb,dim,fobj,fhd)

%%%%-------------------Definitions--------------------------%%
%%
Sun_Pos=zeros(1,dim); % A vector to include the best-so-far Solution, representing the Sun
Sun_Score=inf; % A Scalar variable to include the best-so-far score
Convergence_curve=zeros(1,Tmax);

%%-------------------Controlling parameters--------------------------%%
%%
Tc=3;
M0=0.1;
lambda=15;
%% Step 1: Initialization process
%%---------------Initialization----------------------%%
 % Orbital Eccentricity (e)   
orbital=rand(1,SearchAgents_no); %% Eq.(4)
 %% Orbital Period (T) 
T=abs(randn(1,SearchAgents_no)); %% Eq.(5)
Positions=initialization(SearchAgents_no,dim,ub,lb); % Initialize the positions of planets
t=0; %% Function evaluation counter 
%%
%%---------------------Evaluation-----------------------%%
for i=1:SearchAgents_no
    %% Test suites of CEC-2014, CEC-2017, CEC-2020, and CEC-2022
    PL_Fit(i)=feval(fhd, Positions(i,:)',fobj);
    % Update the best-so-far solution
    if PL_Fit(i)<Sun_Score % Change this to > for maximization problem
       Sun_Score=PL_Fit(i); % Update the best-so-far score
       Sun_Pos=Positions(i,:); % Update te best-so-far solution
    end
end

while t<Tmax %% Termination condition  
 [Order] = sort(PL_Fit);  % Sorting the fitness values of the solutions in current population
 %% The worst Fitness value at function evaluation t
 worstFitness = Order(SearchAgents_no); %% Eq.(11)
 M=M0*(exp(-lambda*(t/Tmax))); %% Eq. (12)
 %% Computer R that represents the Euclidian distance between the best-so-far solution and the ith solution
 for i=1:SearchAgents_no
    R(i)=0;
    for j=1:dim
       R(i)=R(i)+(Sun_Pos(j)-Positions(i,j))^2; %% Eq.(7)
    end
    R(i)=sqrt(R(i));
 end
 %% The mass of the Sun and object i at time t is computed as follows:
 for i=1:SearchAgents_no
    sum=0;
    for k=1:SearchAgents_no
        sum=sum+(PL_Fit(k)-worstFitness);
    end
    MS(i)=rand*(Sun_Score-worstFitness)/(sum); %% Eq.(8)
    m(i)=(PL_Fit(i)-worstFitness)/(sum); %% Eq.(9)
 end
 
 %% Step 2: Defining the gravitational force (F)
 % Computing the attraction force of the Sun and the ith planet according to the universal law of gravitation:
 for i=1:SearchAgents_no
    Rnorm(i)=(R(i)-min(R))/(max(R)-min(R)); %% The normalized R (Eq.(24))
    MSnorm(i)=(MS(i)-min(MS))/(max(MS)-min(MS)); %% The normalized MS
    Mnorm(i)=(m(i)-min(m))/(max(m)-min(m)); %% The normalized m
    Fg(i)=orbital(i)*M*((MSnorm(i)*Mnorm(i))/(Rnorm(i)*Rnorm(i)+eps))+(rand); %% Eq.(6)
 end
 %% a1 represents the semimajor axis of the elliptical orbit of object i at time t,
 for i=1:SearchAgents_no
    a1(i)=rand*(T(i)^2*(M*(MS(i)+m(i))/(4*pi*pi)))^(1/3); %% Eq.(23)
 end
 
 for i=1:SearchAgents_no
    %% a2 is a cyclic controlling parameter that is decreasing gradually from -1 to ?2
    a2=-1+-1*(rem(t,Tmax/Tc)/(Tmax/Tc)); %% Eq.(29)
    %% ? is a linearly decreasing factor from 1 to ?2
    n=(a2-1)*rand+1; %% Eq.(28)
    a=randi(SearchAgents_no);  %% An index of a solution selected at random
    b=randi(SearchAgents_no); %% An index of a solution selected at random
    rd=rand(1,dim); %% A vector generated according to the normal distribution
    r=rand; %% r1 is a random number in [0,1]
    %% A randomly-assigned binary vector 
    U1=rd<r; %% Eq.(21)  
    O_P=Positions(i,:); %% Storing the current position of the ith solution
    %% Step 6: Updating distance with the Sun
    if rand<rand
         %% h is an adaptive factor for controlling the distance between the Sun and the current planet at time t
         h=(1/(exp(n.*randn))); %% Eq.(27)
         %% An verage vector based on three solutions: the Current solution, best-so-far solution, and randomly-selected solution
         Xm=(Positions(b,:)+Sun_Pos+Positions(i,:))/3.0; 
         Positions(i,:)=Positions(i,:).*U1+(Xm+h.*(Xm-Positions(a,:))).*(1-U1); %% Eq.(26)
    else
        
        %% Step 3: Calculating an object’ velocity
        % A flag to opposite or leave the search direction of the current planet
         if rand<0.5 %% Eq.(18)
           f=1;
         else
           f=-1;
         end
         L=(M*(MS(i)+m(i))*abs((2/(R(i)+eps))-(1/(a1(i)+eps))))^(0.5); %% Eq.(15)
         U=rd>rand(1,dim); %% A binary vector
         if Rnorm(i)<0.5 %% Eq.(13)
            M=(rand.*(1-r)+r); %% Eq.(16)
            l=L*M*U; %% Eq.(14)
            Mv=(rand*(1-rd)+rd); %% Eq.(20)
            l1=L.*Mv.*(1-U);%% Eq.(19)
            V(i,:)=l.*(2*rand*Positions(i,:)-Positions(a,:))+l1.*(Positions(b,:)-Positions(a,:))+(1-Rnorm(i))*f*U1.*rand(1,dim).*(ub-lb); %% Eq.(13a)
         else
            U2=rand>rand; %% Eq. (22) 
            V(i,:)=rand.*L.*(Positions(a,:)-Positions(i,:))+(1-Rnorm(i))*f*U2*rand(1,dim).*(rand*ub-lb);  %% Eq.(13b)
         end %% End IF
         
         %% Step 4: Escaping from the local optimum
         % Update the flag f to opposite or leave the search direction of the current planet  
         if rand<0.5 %% Eq.(18)
            f=1;
         else
            f=-1;
         end
         %% Step 5
         Positions(i,:)=((Positions(i,:)+V(i,:).*f)+(Fg(i)+abs(randn))*U.*(Sun_Pos-Positions(i,:))); %% Eq.(25)
    end %% End If
    %% Return the search agents that exceed the search space's bounds
    if rand<rand
       for j=1:size(Positions,2)
          if  Positions(i,j)>ub(j)
              Positions(i,j)=lb(j)+rand*(ub(j)-lb(j));
          elseif  Positions(i,j)<lb(j)
              Positions(i,j)=lb(j)+rand*(ub(j)-lb(j));
          end %% End If
       end   %% End For
    else
       Positions(i,:) = min(max(Positions(i,:),lb),ub);
    end %% End If
    %% Test suites of CEC-2014, CEC-2017, CEC-2020, and CEC-2022
    % Calculate objective function for each search agent
    PL_Fit1=feval(fhd, Positions(i,:)',fobj);
 %  Step 7: Elitism, Eq.(30)
    if PL_Fit1<PL_Fit(i) % Change this to > for maximization problem
       PL_Fit(i)=PL_Fit1; % 
       % Update the best-so-far solution
       if PL_Fit(i)<Sun_Score % Change this to > for maximization problem
           Sun_Score=PL_Fit(i); % Update the best-so-far score
           Sun_Pos=Positions(i,:); % Update te best-so-far solution
       end
    else
       Positions(i,:)=O_P;
    end %% End IF
    t=t+1; %% Increment the current function evaluation
    if t>Tmax %% Checking the termination condition
        break;
    end %% End IF
    Convergence_curve(t)=Sun_Score; %% Set the best-so-far fitness value at function evaluation t in the convergence curve 
 end %% End for i
end %% End while 
Convergence_curve(t-1)=Sun_Score;
end%% End Function


