function val = random_distribution(data_length,varargin)

val = zeros(data_length,1);
if nargin == 1
    for i = 1:data_length
        choice = randi(3,1,1);
        if choice == 1
            a = -0.0015;b=0.0015;
        elseif choice == 2
            a = -0.00175;b = 0.00175;
        else
            a = -0.002;b=0.002;
        end
        val(i,1)= a + (b-a).*rand(1,1);  
    end
end
if nargin == 2
    val = -varargin{1} + 2*varargin{1}.*rand(data_length,1);
end
  
end

