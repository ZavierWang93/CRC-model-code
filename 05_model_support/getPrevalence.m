function prevalence = getPrevalence(yy)

prevalence(:,1) = yy(:,1);
for ii=2:size(yy,2)
    prevalence(:,ii) = (1-sum(prevalence(:,1:ii-1),2)).*yy(:,ii);
end
prevalence(:,end+1) = 1-sum(prevalence,2);

end