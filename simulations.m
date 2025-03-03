function [siz] = simulations(experiment,ts)

% Data
Data1=[3 0; 8 0; 26 0; 76 0; 225 9; 298 17; 258 105; 233 162; 189 176; 128 166; 68 150; 29 85; 14 47; 4 20];
Data=[Data1(1:ts,1); Data1(1:ts,2)];

Sam=2*10^3;                      % desired sample size
Q=2*10^7;                        % maximum number of runs
Out_final=[];                    % Output
sch=0;                           % count of runs (rem 1000)
siz=0;                           % count of accepted runs
q=0;                             % count of runs
eps=0.03;                        % tolerance

while  (q<=Q)&(siz<=Sam)
    % Count runs
    q=q+1;
    if rem(q,1000)==0
        sch=sch+1
    end

    % Random numbers generator
    rng('shuffle','twister')
    seed=rng;

    % Sample from priors
    th1=0.008*rand;
    th2=0.5+3*rand;
    th3=1.6*rand;
    th4=1.1*rand;

    % Simulate the stochastic process
    [Output] = Gil3([th1 th2 th3 th4]);

    % Calculate the statistics
    S=[Output(29:29+ts-1);Output(43:43+ts-1)];
    S=(S-Data).^2;
    S=sum(S);
    S_new=S/sum(Data.^2);

    % Accept or reject this run and update the sample size
    if S_new<=eps
        Out_final=[Out_final [th1; th2; th3; th4; Output]];
        siz=siz+1;
        siz
    end
end

% Save the results
writematrix(Out_final,['Out_final' num2str(experiment) num2str(ts) '.csv'])

end