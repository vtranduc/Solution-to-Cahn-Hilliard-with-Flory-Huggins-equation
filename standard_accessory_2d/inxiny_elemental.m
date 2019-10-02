function [inx,iny]=inxiny_elemental(element,ney)
iny=mod(element,ney);
if iny==0
    iny=ney;
end
inx=(element-iny)/ney+1;
end