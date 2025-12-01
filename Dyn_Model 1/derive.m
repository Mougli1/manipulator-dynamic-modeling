function df = derive(t,f)

df=[f(1,2)-f(1,1),(f(1,3:end)-f(1,1:end-2))/2,f(1,end)-f(1,end-1)]/(t(1,2)-t(1,1));

end 