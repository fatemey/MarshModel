function printjobs(n)

for i = 1 : n
    
    fileID = fopen(['job',num2str(i)],'w');
    fprintf(fileID,'#!/bin/bash\n\n');
    fprintf(fileID,'#SBATCH -p scavenger\n');
    fprintf(fileID,'#SBATCH --job-name=marsh%d\n',i);
    fprintf(fileID,'#SBATCH --output=out%d.txt\n',i);
    fprintf(fileID,'#SBATCH --cpus-per-task=32\n');
    fprintf(fileID,'#SBATCH --mem=20000\n');
    fprintf(fileID,'module load Matlab/R2016a\n');
    fprintf(fileID,'matlab -nodisplay -nodesktop -singleCompThread -r "i=%d; n=3^7/%d; Run((i-1)*n+1, i*n)"\n',i,n);
    fclose(fileID);
    
end


fileID = fopen('jobs.sh','w');
for i = 1 : n
    fprintf(fileID,'sbatch job%d\n',i);
end
fclose(fileID);

end