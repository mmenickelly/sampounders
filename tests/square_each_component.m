function [F,J] = square_each_component(func,x,Set)
    
    [Ftemp,Jtemp] = func(x,Set);
    F = 0.5*Ftemp.^2;
    
    n = length(x);
    J = Jtemp.*repmat(Ftemp,n,1);
end
    
    