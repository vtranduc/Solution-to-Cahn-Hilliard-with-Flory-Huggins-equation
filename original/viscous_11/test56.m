function test56

hey=foo(1)

end

function sol=foo(arg)

if arg<5
    take=arg+1;
    sol=foo(take);
end

end