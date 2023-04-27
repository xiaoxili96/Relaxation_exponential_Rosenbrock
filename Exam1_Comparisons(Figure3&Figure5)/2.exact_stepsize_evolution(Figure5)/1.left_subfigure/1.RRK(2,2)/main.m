stepsize=[1/2 1/4 1/8];

for k=1:size(stepsize,2)
    tau=stepsize(1,k);
    [exact_stepsize,tmesh]=rrk2s2(tau);
    eval([['exact_stepsize_rrk2s2_' num2str(1/tau)] '=exact_stepsize;']);
    eval([['tmesh_rrk2s2_' num2str(1/tau)] '=tmesh;']);
end

clear stepsize k exact_stepsize tmesh tau
save('data_rrk2s2.mat');
