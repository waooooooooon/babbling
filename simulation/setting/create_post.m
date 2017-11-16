%initialization
Position_xy = zeros(1000,2);
Post_position = zeros(1000,100);
Post_ID = zeros(1000,100);

nrows = 100;
ncols = 10;
NeuronID_Position = reshape(randperm(ncols*nrows), [nrows ncols]);

r = 1;


%create Position information
k = 1;
for i = 1:10
    for j = 1:100
       
        Position_xy(k,:) = [i,j];
        k = k + 1;
    end
end

%caliculate the distance and nearly neuron ID
[idx, dist] = rangesearch(Position_xy,Position_xy,r*sqrt(2));

%caliculate Post neuron 
for i = 1: 1000
    size_idx = size(idx{i});
    Post_position(i,1:size_idx(1,2)) = idx{i};
end
%delete number of themselve
Post_position(:,1) = 0;


for i = 1:1000
size_neuron = 2:max(find(Post(find(NeuronID_Position==i),:)~=0));
Post_ID(i,size_neuron) = NeuronID_Position(Post_position(find(NeuronID_Position==i),size_neuron)); %this sentense create new post matrix
end


% neuron caliculation


Position;
