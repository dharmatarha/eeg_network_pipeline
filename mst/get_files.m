%function [datastruct]=get_files(dirname)

cd(dirname)
files=dir('*.txt');

files=rmfield(files, {'date','bytes','isdir','datenum'});


for x=1:length({files.name})
files(x).PLI=textread(files(x).name);
files(x).RND_PLI=randomize_PLI_matrix(files(x).PLI);
end

mean_PLI=mean(reshape([files.PLI],[size(files(1).PLI),length({files.name})]),3);
mean_RND_PLI=mean(reshape([files.RND_PLI],[size(files(1).RND_PLI),length({files.name})]),3);