function n=determine_nIterations(folder)
n=0;
while 1==1
    n=n+1;
    if ~exist([folder 'concentration/iteration_' num2str(n)],'file')
        break
    end
end
n=n-1;
end