close all; clear all; clc;
%% Problem description
disp('=============================================================')
disp('Genetic algorithms: the Knapsack Problem')
% By: Darwineswaran Raja Lingam
disp('=============================================================')
disp('=============================================================================')
disp('Problem: Bob went into a room with bag of 25kg max. capacity. He has been    ')
disp('         presented 10 items with various weights and values. He needs to ')
disp('         collect the items with as large value as possible without exceeding ')
disp('         the mass capacity. So what is the best combination of items for Bob?')
disp('=============================================================================')
%% List of items
item_name = {'Belt';'Ball';'Apple';'Rat';'Dog';'Cat';'Box';'Pen';'Door';'Book'};
WEIGHT = [3; 10; 2; 3; 18; 5; 13; 6; 10; 5;];
Value = [2; 5; 6; 10; 18; 13; 10; 8; 14; 20;];
item_list = table(WEIGHT, Value, 'RowNames', item_name);
disp(item_list)
%% Genetic algorithm
% Problem initialization
n = 10;   % items
Weight = [3 10 2 3 18 5 13 6 10 5];    % Weight
Value = [2 5 6 10 18 13 10 8 14 20];    % Value
Capacity = 25;   % Capacity
disp('=============================================================================')
disp('By using brute force method, maximum value and the best combination obtained')
[bF_maxValue, bF_bestCombination] = bruteForce(Weight, Value, n, Capacity) % bruteForce value (control)
disp('=============================================================================')
popsize = 50;   % Population size

popbest = zeros(popsize, n);   % Optimal population storage object

pop = initpop(popsize,n); % Initial population

gen = 1000;   % Number of generations

pc = 0.9;   % Crossover rate
pm = 0.09;  % Mutation rate

% Record the most suitable individuals of the population for each iteration
% (graph)
y = zeros(1,popsize);   % Record maximum value
g = zeros(1,popsize);   % Record maximum weight
n = zeros(1,popsize);   % Record location

% Fitness calculation
for i = 1:gen
    
    [fitvalue] = calobjvalue(pop,Value,Weight,Capacity); %Calculate the objective function
                                        %Correspondence measure of environmental fitness
    
    [bestweight,bestvalue, bestpop] = best(pop,fitvalue,Weight); 
    %Calculate optimal individual weight, value, and location
    disp("Generation = " + i + ", Fitness value = " + bestvalue)
    y(i) = max(bestvalue);    %Record maximum value
    
    g(i) = max(bestweight);   %Record maximum weight
    
    n(i) = i;                 %Record location
    
    popbest(i,:) = bestpop;    %Record the best individual      
    
    [newpop] = selection(pop,fitvalue);  %Perform selection calculations
                                       %Corresponding to find individuals who have the right to reproduce                                
    [newpop1] = crossover(newpop,pc);    %Cross calculation
                                       %Corresponding gene exchange                             
    [newpop2] = mutation(newpop1,pm);     %Perform variation calculations
                                        %Corresponding gene mutation
    pop = newpop2;   %Iterative repetition
end
pop;

% Graph visualization
i = 1:gen;
plot(y(i),'-b')
xlabel('Generation')
ylabel('Fitness value');
title('The Genetic Algorithm optimization result');
legend('Mean');
grid on

xlim([0 gen])
ylim([25 70]) % limit based on item value

[z, index] = max(y);
po = n(index);         %Optimal algebraic position
GA_maxWeight = g(index)         %Optimal weight
disp('=============================================================================')
disp('By using genetic algorithm, maximum value and the best combination obtained')
GA_maxValue = z                %Best value
GA_bestCombination = popbest(po,:)    %Optimal solution
disp('=============================================================================')
%% Crossover function
function [newpop1] = crossover(newpop,pc)

[popsize,individual] = size(newpop);
newpop1 = zeros(popsize,individual);

for i = 1:2:popsize-1
    ps = rand;
    
    if ps < pc
        cpoint = round(rand*individual);
        
        newpop1(i,:) = [newpop(i,1:cpoint),newpop(i+1,cpoint+1:individual)];
        
        newpop1(i+1,:) = [newpop(i+1,1:cpoint),newpop(i,cpoint+1:individual)];
    else
        newpop1(i,:) = newpop(i,:);
        newpop1(i+1,:) = newpop(i+1,:);
    end
end
end
%% Mutation function
function [newpop2] = mutation(newpop1,pm)

[popsize,individual] = size(newpop1);

for i = 1:popsize
    ps = rand;
    if ps < pm
        mpoint = round(rand*individual);
        if mpoint <= 0
            mpoint = 1;
        end
        if newpop1(i,mpoint) == 0
            newpop1(i,mpoint) = 1;
        else
            newpop1(i,mpoint) = 0;
        end
    else
        
    end
end
newpop2 = newpop1;
end
%% Selection Function
function [newpop] = selection(pop,fitvalue)

    totalfit = sum(fitvalue);
    
    pfitvalue = fitvalue/totalfit; 
  
    mfitvalue = cumsum(pfitvalue); 
    
    [popsize,~] = size(pop); 
    
    ms = sort(rand(popsize,1)); 
    
    fitin = 1; 
    newin = 1;
    newpop = zeros(size(pop));
    
    while newin <= popsize
        if mfitvalue(fitin) > ms(newin)
            newpop(newin,:) = pop(fitin,:); 
            newin = newin+1; 
        else 
        fitin = fitin + 1; 
        end 
    end
end
%% Initialize population
function pop = initpop(popsize,n)
pop = round(rand(popsize,n));
end
%% Fitness function calculation
% Calculate the fitness of each individual, and use the penalty function to calculate the objective function
% Using the penalty function method to calculate the objective function and its fitness
function [fitvalue] = calobjvalue(pop,V,W,CW)

[popsize,individual] = size(pop);

sumV = zeros(1,1000);
sumW = zeros(1,1000);
fitvalue = zeros(1,1000);

for i = 1:popsize
    for j = 1:individual
          sumV(i) = sumV(i)+V(j) * pop(i,j);   % Calculate the individual value and as fitness
          
          sumW(i) = sumW(i)+W(j) * pop(i,j);   % Calculate the weight of the individual as a penalty function
                                           % Not exceeding the boundary is itself, exceeding the boundary infinity
    end
    if sumW(i) > CW
          sumV(i) = 0;
    else 
       sumV(i) = sumV(i);
    end
  fitvalue(i) = sumV(i);
end
end
function [bestweight,bestvalue,bestpop] = best(newpop,fitvalue,w)
[popsize,individual] = size(newpop);
bestvalue = fitvalue(1);
for i = 2:popsize
    if fitvalue(i) > bestvalue
        bestvalue = fitvalue(i); 
    end
end
[~,index] = max(fitvalue);

bestweight = 0;
bestpop = zeros(1,individual);
i = index;
for j = 1:individual
    bestweight = w(j) * newpop(i,j) + bestweight;
    bestpop(1,j) = newpop (i,j);
end
end
%% Brute force function
function [maxValue, maxComb] = bruteForce(Weight, Profit, n, capacity)

items = 1:n;

maxValue = 0;
maxComb = [];

for i = items
    combination = nchoosek(items,i);
    
    for row = 1:size(combination)
        tmpV = 0;
        tmpW = 0;
        t = combination(row,:);
        for elem = t
            tmpW = tmpW + Weight(elem);
            tmpV = tmpV + Profit(elem);            
        end
        if tmpW <= capacity && tmpV>maxValue
            maxValue = tmpV;
            maxComb = Weight(t);
        end
    end
end

maxValue;
maxComb;
end
