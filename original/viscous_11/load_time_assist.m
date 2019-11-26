function sol=load_time_assist(keys,load,load_more,folder)
sol=zeros(1,load_more);
for i=1:1:load
    sol(i)=dlmread([folder 'time/iteration_' num2str(keys(i))]);
end
end