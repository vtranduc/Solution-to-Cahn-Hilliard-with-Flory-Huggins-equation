function solution=get_elements_3d(xth,yth,zth,nex,ney,nx,ny,nz,...
    positional_element)
%positional_element is an optional variable. If it is not specified,
%all adjacent elements will be computed. It will also test for existence
%of adjacent elements.
%If positional_element is specified, it WILL NOT check if the
%positional_element actually exists. So, you have to make sure
%it exists.
if nargin==8
    solution=zeros(1,8);
    z1=nex*ney*(zth-1);
    z2=nex*ney*(zth-2);
    y1=nex*(yth-1);
    y2=nex*(yth-2);
    x1=(xth-1);
    if xth==1
        if yth==1
            if zth==1
                solution(8)=z1+y1+xth;
            elseif zth==nz
                solution(4)=z2+y1+xth;
            else
                solution(4)=z2+y1+xth;
                solution(8)=z1+y1+xth;
            end
        elseif yth==ny
            if zth==1
                solution(6)=z1+y2+xth;
            elseif zth==nz
                solution(2)=z2+y2+xth;
            else
                solution(2)=z2+y2+xth;
                solution(6)=z1+y2+xth;
            end
        else
            if zth==1
                solution(6)=z1+y2+xth;
                solution(8)=z1+y1+xth;
            elseif zth==nz
                solution(2)=z2+y2+xth;
                solution(4)=z2+y1+xth;
            else
                solution(2)=z2+y2+xth;
                solution(4)=z2+y1+xth;
                solution(6)=z1+y2+xth;
                solution(8)=z1+y1+xth;
            end
        end
    elseif xth==nx
        if yth==1
            if zth==1
                solution(7)=z1+y1+x1;
            elseif zth==nz
                solution(3)=z2+y1+x1;
            else
                solution(3)=z2+y1+x1;
                solution(7)=z1+y1+x1;
            end
        elseif yth==ny
            if zth==1
                solution(5)=z1+y2+x1;
            elseif zth==nz
                solution(1)=z2+y2+x1;
            else
                solution(1)=z2+y2+x1;
                solution(5)=z1+y2+x1;
            end
        else
            if zth==1
                solution(5)=z1+y2+x1;
                solution(7)=z1+y1+x1;
            elseif zth==nz
                solution(1)=z2+y2+x1;
                solution(3)=z2+y1+x1;
            else
                solution(1)=z2+y2+x1;
                solution(3)=z2+y1+x1;
                solution(5)=z1+y2+x1;
                solution(7)=z1+y1+x1;
            end
        end
    else
        if yth==1
            if zth==1
                solution(7)=z1+y1+x1;
                solution(8)=z1+y1+xth;
            elseif zth==nz
                solution(3)=z2+y1+x1;
                solution(4)=z2+y1+xth;
            else
                solution(3)=z2+y1+x1;
                solution(4)=z2+y1+xth;
                solution(7)=z1+y1+x1;
                solution(8)=z1+y1+xth;
            end
        elseif yth==ny
            if zth==1
                solution(5)=z1+y2+x1;
                solution(6)=z1+y2+xth;
            elseif zth==nz
                solution(1)=z2+y2+x1;
                solution(2)=z2+y2+xth;
            else
                solution(1)=z2+y2+x1;
                solution(2)=z2+y2+xth;
                solution(5)=z1+y2+x1;
                solution(6)=z1+y2+xth;
            end
        else
            if zth==1
                solution(5)=z1+y2+x1;
                solution(6)=z1+y2+xth;
                solution(7)=z1+y1+x1;
                solution(8)=z1+y1+xth;
            elseif zth==nz
                solution(1)=z2+y2+x1;
                solution(2)=z2+y2+xth;
                solution(3)=z2+y1+x1;
                solution(4)=z2+y1+xth;
            else
                solution(1)=z2+y2+x1;
                solution(2)=z2+y2+xth;
                solution(3)=z2+y1+x1;
                solution(4)=z2+y1+xth;
                solution(5)=z1+y2+x1;
                solution(6)=z1+y2+xth;
                solution(7)=z1+y1+x1;
                solution(8)=z1+y1+xth;
            end
        end
    end
else
    if positional_element==1
        solution=nex*ney*(zth-2)+nex*(yth-2)+(xth-1);
    elseif positional_element==2
        solution=nex*ney*(zth-2)+nex*(yth-2)+xth;
    elseif positional_element==3
        solution=nex*ney*(zth-2)+nex*(yth-1)+(xth-1);
    elseif positional_element==4
        solution=nex*ney*(zth-2)+nex*(yth-1)+xth;
    elseif positional_element==5
        solution=nex*ney*(zth-1)+nex*(yth-2)+(xth-1);
    elseif positional_element==6
        solution=nex*ney*(zth-1)+nex*(yth-2)+xth;
    elseif positional_element==7
        solution=nex*ney*(zth-1)+nex*(yth-1)+(xth-1);
    elseif positional_element==8
        solution=nex*ney*(zth-1)+nex*(yth-1)+xth;
    end
end
end

% solution(1)=z2+y2+x1;
% solution(2)=z2+y2+xth;
% solution(3)=z2+y1+x1;
% solution(4)=z2+y1+xth;
% solution(5)=z1+y2+x1;
% solution(6)=z1+y2+xth;
% solution(7)=z1+y1+x1;
% solution(8)=z1+y1+xth;