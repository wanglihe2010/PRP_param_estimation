function demand2 = pulse_aggreation(time, int, dur, sim_s, T)
n_new = length(find(time<=sim_s));
%demand = zeros(sim_s/T,1);
demand2 = zeros(sim_s,1);
for i = 1:n_new
    demand2(time(i):min(time(i) + dur(i)-1,sim_s)) = demand2(time(i):min(time(i) + dur(i)-1,sim_s)) + int(i);
end
if T >1  
    demand2 = sum(reshape(demand2,T, length(demand2)/T))';
end

%if T == 1
   % for i = 1:n_new
   %     demand(time(i):min(time(i) + dur(i)-1,sim_s/T)) = demand(time(i):min(time(i) + dur(i)-1,sim_s/T)) + int(i);
   % end
%else 

 %   index1 = ceil(time/T);
  %  index2 = ceil((time+ dur)/T);
   % %index2 = min(ceil((time+dur)/T),sim_s/T);
    %if max(index2) > sim_s/T
     %   demand = zeros(max(index2),1);
   % end
    %for i = 1:n_new


        % The follow puts the demands into the appropriate 1-hour "bins"
    %    if index1(i)==index2(i)  % if completely within one 1-hour time period
     %       demand(index1(i),1) = demand(index1(i),1) + dur(i)*int(i); 
      % else
       %     for ii = index1(i):index2(i)  % splits the demands between the appropriate hourly bins
        %        if ii == index1(i)
         %           if ii == 0
          %              ii=1;
           %         end
            %        demand(ii,1) = demand(ii,1) + (ceil(time(i)/T)-(time(i)-1)/T)*int(i)*T;
             %   elseif ii == index2(i)
              %      if ii >length(demand)
               %         ii = 99361;
                %    end
                 %   demand(ii,1) = demand(ii,1) + ((time(i)+dur(i)-1)/T-floor((time(i)+dur(i)-1)/T))*int(i)*T;
                %else
                 %   demand(ii,1) = demand(ii,1) + int(i)*T;
                %end
           %end
        %end

    %end
    %demand = demand(1:sim_s/T,:);
%end