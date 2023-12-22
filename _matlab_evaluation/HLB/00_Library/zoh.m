function Vq = zoh(t,V,tq)
j = 1;
for i=1:length(tq)
    if (tq(i) > t(j))
        j = j + 1;
    end
    Vq(i) = V(j);
end
end