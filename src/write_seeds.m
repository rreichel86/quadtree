function write_seeds(Filename,seeds)

num_seeds = size(seeds,1);
dim = size(seeds,2);

if dim == 2
    fileID = fopen(Filename,'w');
    fprintf(fileID,'%6d\n',num_seeds);
    for i = 1: num_seeds
        fprintf(fileID,'%32.16f %32.16f\n',seeds(i,:));
    end 
    fclose(fileID);

elseif dim == 3
    fileID = fopen(Filename,'w');
    fprintf(fileID,'%6d\n',num_seeds);
    for i = 1: num_seeds
        fprintf(fileID,'%32.16f %32.16f %32.16f\n',seeds(i,:));
    end 
    fclose(fileID);
end

end
