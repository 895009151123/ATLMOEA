function Offspring = OperatorDE1(Parent1,Parent2,Parent3,Lx1,t,Parameter)
    %% Parameter setting
%     if nargin > 3
%         [CR,F,proM,disM] = deal(Parameter{:});
%     else
        [CR,F,proM,disM] = deal(1,0.5,1,20);
%     end
    if isa(Parent1(1),'SOLUTION')
        calObj  = true;
        Parent1 = Parent1.decs;
        Parent2 = Parent2.decs;
        Parent3 = Parent3.decs;
    else
        calObj = false;
    end
    [N,D]   = size(Parent1);
    Problem = PROBLEM.Current();

    %% Differental evolution
    Site = rand(N,D) < CR;
    Offspring       = Parent1;
    Offspring(Site) = Offspring(Site) + t*(Parent2(Site)-Parent3(Site))+ (1-t)*(Parent1(Site)-Lx1(Site));

    %% Polynomial mutation
    Lower = repmat(Problem.lower,N,1);
    Upper = repmat(Problem.upper,N,1);
    Site  = rand(N,D) < proM/D;
    mu    = rand(N,D);
    temp  = Site & mu<=0.5;
    Offspring       = min(max(Offspring,Lower),Upper);
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                      (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5; 
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                      (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
    if calObj
        Offspring = SOLUTION(Offspring);
    end
end