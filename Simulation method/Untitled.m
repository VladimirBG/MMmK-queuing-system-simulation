clc, clear, close
seed = 456;
rng('shuffle');
Creq = 0;
    ttemp = 0;
     for J = 1 : 50
        ttemp = (-1/0.05*log(rand));
        treq(J) = ttemp;
    end
    lambdahat = poissfit(treq);
    
    display(lambdahat);
    
