function printjobs(m,n) % m: number of total samples, n: number of groups

for i = 1 : n
    fileID = fopen(['job',num2str(i)],'w');
    fprintf(fileID,'#!/bin/bash\n\n');
    fprintf(fileID,'#SBATCH -p scavenger\n');
    fprintf(fileID,'#SBATCH --job-name=marsh%d\n',i);
    fprintf(fileID,'#SBATCH --output=out%d.txt\n',i);
    fprintf(fileID,'#SBATCH --cpus-per-task=31\n');
    fprintf(fileID,'#SBATCH --mem=20000\n');
    fprintf(fileID,'module load Matlab/R2016a\n');
    fprintf(fileID,'matlab -nodisplay -nodesktop -singleCompThread -r "i=%d; n=floor(%d/%d); Run((i-1)*n+1, i*n)"\n',i,m,n);
    fclose(fileID);
end

if rem(m,n) ~= 0
    fileID = fopen(['job',num2str(n)],'w');
    fprintf(fileID,'#!/bin/bash\n\n');
    fprintf(fileID,'#SBATCH -p scavenger\n');
    fprintf(fileID,'#SBATCH --job-name=marsh%d\n',n);
    fprintf(fileID,'#SBATCH --output=out%d.txt\n',n);
    fprintf(fileID,'#SBATCH --cpus-per-task=31\n');
    fprintf(fileID,'#SBATCH --mem=20000\n');
    fprintf(fileID,'module load Matlab/R2016a\n');
    fprintf(fileID,'matlab -nodisplay -nodesktop -singleCompThread -r "i=%d-rem(%d,%d)-floor(%d/%d)+1; n=%d; Run(i,n)"\n',m,m,n,m,n,m);
    fclose(fileID);
end

fileID = fopen('jobs.sh','w');
for i = 1 : n
    fprintf(fileID,'read -st 2; sbatch job%d\n',i);
end
fclose(fileID);

end