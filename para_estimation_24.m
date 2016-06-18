function [para_est] = para_estimation_24(demand1,demand2,t1,t2)
%%%output unit: lambda(/s), alpha, beta(l or g/s), tau(s)



demand1_reshape = reshape(demand1,3600/t1*24,length(demand1)*t1/3600/24)';
demand1_reshape = reshape(demand1_reshape,length(demand1)/24,24);
demand2_reshape = reshape(demand2,3600/t2*24,length(demand2)*t2/3600/24)';
demand2_reshape = reshape(demand2_reshape,length(demand2)/24,24);
for j = 1:24 
    para_est(j,:) = para_estimation(demand1_reshape(:,j),demand2_reshape(:,j),t1,t2);
end

end
