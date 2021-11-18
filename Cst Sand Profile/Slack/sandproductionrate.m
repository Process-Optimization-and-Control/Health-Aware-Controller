function [sandArray,sandArrayNoise,stepSandArray] = sandproductionrate(initial_mdot,days,type,k)

if (type == 'exp')
sandArray = [initial_mdot];

    for n = 1:days
     %sandArray = [sandArray initial_mdot*(0.5+exp(k*n))];
     sandArray = [sandArray initial_mdot*(exp(k*n))];
    end 

    stepSandArray = [];
    step = 1;

    for n = 1:days
       stepSandArray = [stepSandArray sandArray(step)];
       if mod(n,30) == 0
           step = n;

       end

    end
end

if (type == 'log')
sandArray = (10*initial_mdot)/(1+exp(-k*(0-days/2)));

    for n = 1:days
     newVal = (10*initial_mdot)/(1+exp(-k*(n-days/2)));
     sandArray = [sandArray newVal];
    end 

    stepSandArray = [];
    step = 1;

    for n = 1:days
       stepSandArray = [stepSandArray sandArray(step)];
       if mod(n,30) == 0
           step = n;

       end

    end
end

if (type == 'cst')
    sandArray = initial_mdot*ones(1,days+1);
end

sandArrayNoise = sandArray + 0.1 * rand(1, length(sandArray));
end

