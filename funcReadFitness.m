function [fitness,fitnessBaseline,NoG_actual] = funcReadFitness(caseName,runName,NoNests,Ma)

fid = fopen([caseName,runName,'Fitness',num2str(Ma,1),'.txt']);

if fid == -1
    error('ERROR(funcReadFitness): Fitness file not found.');
end

% Read a number of columns equal to number of nests.
readString = '%d';
for i = 1:NoNests
    readString = strcat(readString,'%f');
end
fitnessRead = textscan(fid,readString,'HeaderLines',1);
NoNests = length(fitnessRead)-1;
NoG_actual = length(fitnessRead{1});
fclose(fid);

% Format as matrix.
fitness = zeros(NoG_actual, NoNests);
for i = 2:NoNests+1
    fitness(:,i-1) = [fitnessRead{i}];
end

% Read baseline fitness.
fid = fopen([caseName,runName,'Fitness_0.txt']);
fitnessBaseline = fscanf(fid,'%f');
fitnessBaseline = fitnessBaseline(2);

% Print fitness.
fprintf('\nOPTIMISER RESULTS:\n');
fprintf('       BASELINE FITNESS: %f\n',fitnessBaseline);
fitnessOptimal = fitness(NoG_actual,1);
fprintf('      OPTIMISED FITNESS: %f\n',fitnessOptimal);
fitnessIncrease = (1-fitnessBaseline/fitnessOptimal)*100;
fprintf('    INCREASE IN FITNESS: %.2f%%\n\n',fitnessIncrease);

fclose(fid);

end