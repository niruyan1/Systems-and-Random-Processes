function [x,n]=unitstep(n0,n1,n2)
n=n1:n2;
x=ones(size(n));
counter=1;
for b=n1:n2
    if b<n0
        x(counter)=0;
    end
    counter=counter+1;
end
end