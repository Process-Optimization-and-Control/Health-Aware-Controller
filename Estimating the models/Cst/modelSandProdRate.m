%% Model sand production rate with interpolation

function modelSandArray = modelSandProdRate(sandArray,dt)

modelSandArray = [];

days = length(sandArray);

modelSandArray(1) = sandArray(1);
deltaSand = 0;

for i = 2:days
    
    if i <= dt
       modelSandArray(i) = sandArray(i);
    end 
    
    modelSandArray(i) = modelSandArray(i-1);
    
    if mod(i,dt) == 0
        deltaSand = (sandArray(i) - sandArray(i-dt+1))/dt;
    end
    
    modelSandArray(i) = modelSandArray(i-1) + deltaSand;
    
    if mod(i,dt) == 0
        modelSandArray(i) = sandArray(i);
    end
    
end



figure(1)
plot(1:500,modelSandArray(1:500),'b')
hold on
plot(1:500,sandArray(1:500),'r')
legend({'Model Sand Prod Rate','Sample Sand Prod Rate'},'Location','northwest')
title('Sand production rates, modelled, sampled and true')

end
