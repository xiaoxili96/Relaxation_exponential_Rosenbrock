stepsize=[1/2 1/4 1/8];

for k=1:size(stepsize,2)
    tau=stepsize(1,k);
    [exact_stepsize,tmesh]=rexprb3s2(tau);
    eval([['exact_stepsize_rexprb3s2_' num2str(1/tau)] '=exact_stepsize;']);
    eval([['tmesh_rexprb3s2_' num2str(1/tau)] '=tmesh;']);
end

clear stepsize k exact_stepsize tmesh tau
save('data_rexprb3s2.mat');
