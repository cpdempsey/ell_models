function convolved_matrix = convolve_matrix_by_tau(dt,tau,input_matrix)
%takes dt,tau,input_matix, where each row of input_matrix is a spike
%train/time series. Convolves rows of input_matrix with an exponential
%kernel with time constant tau; attenuates result to have same dimensions
%as original

tran=0:dt:200;
kernel=(1/tau)*exp(-tran/tau);

input_padded=[input_matrix(:,floor(end/2):end) input_matrix];

convolved_matrix=zeros(size(input_padded,1),size(input_padded,2)+length(tran)-1);
for i=1:size(input_padded,1)
    convolved_matrix(i,:)=conv(input_padded(i,:),kernel)*dt;
end

convolved_matrix = convolved_matrix(:,...
            floor(size(input_matrix,2)/2+1)+1:(floor(size(input_matrix,2)/2+1)...
                                    + size(input_matrix,2)));
                                
 