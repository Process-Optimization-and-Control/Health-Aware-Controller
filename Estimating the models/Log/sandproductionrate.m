function [sandArray,sandArrayNoise,stepSandArray] = sandproductionrate(initial_mdot,days,type,k,samplingrate)

sandArray = [initial_mdot];


if (type == 'exp')
    for n = 1:days
     sandArray = [sandArray initial_mdot*(0.5+exp(k*n))];
    end 
    
    sandArray = sandArray(2:501);
    
    stepSandArray = [];
    step = 1;

    for n = 1:days
       stepSandArray = [stepSandArray sandArray(step)];
       if mod(n,samplingrate) == 0
           step = n;

       end

    end
end

if (type == 'log')
    for n = 1:days
     newVal = (10*initial_mdot)/(1+exp(-k*(n-days/2)));
     sandArray = [sandArray newVal];
    end 

    sandArray = sandArray(2:501);
    
    stepSandArray = [];
    step = 1;
    
    for n = 1:days
       stepSandArray = [stepSandArray sandArray(step)];
       if mod(n,samplingrate) == 0
           step = n;

       end

    end
end

sandArrayNoise = sandArray + 0.01 * rand(1, length(sandArray));

end
