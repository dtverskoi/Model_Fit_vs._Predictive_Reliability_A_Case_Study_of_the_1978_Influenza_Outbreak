function [Output] = Gil3(theta)

% Parameters
Tmax=30;                             % total time

% Random numbers generator
rng('shuffle','twister')
seed=rng;

% Parameters
S=[-1 0 0 0; 1 -1 0 0; 0 1 -1 0; 0 0 1 -1; 0 0 0 1];

% Outcome
Out=[];

% Initialization
t=0;                                 % time
x=[762; 1; 0; 0; 0];                 % number of molecules
Out=[Out [t; x]];                    % Output
delay=[];

while t<Tmax
    h=NaN(4,1);                      % hazard functions
    h(1)=theta(1)*x(1)*x(2);
    h(2)=theta(2)*x(2);
    h(3)=theta(3)*(x(3)-length(delay));  % some B students are reserved for delayed reactions
    h(4)=theta(4)*x(4);
    h0=sum(h);                       % combined hazard
    if h0~=0
        dt=exprnd(1/h0);                % dt
        if ~isempty(delay) && delay(1) < t+dt     % check if there is a delayd reaction in the set [t,t+dt)
            t=delay(1);                           % update t to be the time of the next delayed reaction
            x=x+S(:,3);
            delay(1) = [];                        % remove the processed reaction
            Out=[Out [t; x]];                     % Output
        else
            t=t+dt;                               % update time
            j=randsample(4,1,true,h/h0);          % determine the reaction index
            if j~=3
                x=x+S(:,j);                       % update the number of molecules
            else
                delay=[delay t+1.5];
            end
            Out=[Out [t; x]];                     % Output
        end
    else
        t=Tmax;
        Out=[Out [t; x]];                    % Output
    end
end

Ou=NaN(6,14);
for k=1:14
    Aux_O=Out(1,:);
    J1=find(Aux_O<=k);
    z1=J1(end);
    Ou(:,k)=Out(:,z1);
end

Z=763-Out(2,end);

Output=[Ou(2,:)'; Ou(3,:)'; Ou(4,:)'; Ou(5,:)'; Z];

end