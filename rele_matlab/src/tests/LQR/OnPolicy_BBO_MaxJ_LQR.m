%% LQR gradient test
addpath(genpath('../../Statistics'));
addpath(genpath('../../iodata'));
addpath('../');

%clear old data
clear all;
clc;
close all;

algorithms{1} = 'pgpe';
algorithms{2} = 'nes';
algorithms{3} = 'rwr';
algorithms{4} = 'enes';
algorithms{5} = 'reps';


distributions{1} = 'gauss';
distributions{2} = 'chol';
distributions{3} = 'diag';
distributions{4} = 'log';

nbEpisodes = 50;
nbUpdates  = 200;
stepLength = 0.01;

domain = 'lqr';

%prog = ['/home/matteo/Projects/github/ReLe/rele-build/',domain,'_BBO'];
prog = ['/home/dave/ReLe/devel/lib/rele/',domain,'_BBO'];

args = [num2str(nbUpdates), ' ', num2str(nbEpisodes), ...
    ' ', num2str(stepLength)];

% figure(1);
% hold on;
nbtests = length(distributions)*2+2;
J = zeros(nbUpdates,nbtests);
count = 1;
for i = 1 : 3
    
    for k = 1 : length(distributions)
        cmd = [prog, ' ', algorithms{i}, ' ', distributions{k}, ' ', args];

        disp('=======================================================================');
        disp(['### Running algorithm: ',algorithms{i}, ' (',  distributions{k},') ###'])
        disp('=======================================================================');
        status = system(cmd);
        
        %% show results
        disp('Reading agent data...')
        csv = csvread(['/tmp/ReLe/',domain,'/BBO/',domain,'_',algorithms{i},'_', distributions{k},'_agentData.log']);
        
        disp('Organizing data...')
        
        index = 1;
        ep = 1;
        
        clear data
        
        if strcmp(algorithms{i},'nes') || strcmp(algorithms{i},'enes')
            while(index < size(csv, 1))
                [data(ep), index] = ReadNESStatistics(csv, index);
                ep = ep + 1;
            end
        elseif strcmp(algorithms{i},'pgpe')
            while(index < size(csv, 1))
                [data(ep), index] = ReadPGPEStatistics(csv, index);
                ep = ep + 1;
            end
        elseif strcmp(algorithms{i},'rwr')
            while(index < size(csv, 1))
                [data(ep), index] = ReadEMStatistics(csv, index);
                ep = ep + 1;
            end
            
        end
        
        clearvars csv
        
        
        J_history = [];
        for o = 1:length(data)
            J(o,count) = mean([data(o).policies.J]);
            %             J_history = [J_history, data(k).J'];
        end
        testname{count} = [algorithms{i}, ' ', distributions{k}];
        count = count + 1;
    end
   
end
for i = 4 : 5
    cmd = [prog, ' ', algorithms{i}, ' ', args];
    status = system(cmd);
    
    if (status==0)
        %% show results
        disp('Reading agent data...')
        csv = csvread(['/tmp/ReLe/',domain,'/BBO/',domain,'_',algorithms{i},'_agentData.log']);
        
        disp('Organizing data...')
        
        index = 1;
        ep = 1;
        
        clear data
        
        if strcmp(algorithms{i},'nes') || strcmp(algorithms{i},'enes')
            while(index < size(csv, 1))
                [data(ep), index] = ReadNESStatistics(csv, index);
                ep = ep + 1;
            end
        elseif strcmp(algorithms{i},'reps')
            while(index < size(csv, 1))
                [data(ep), index] = ReadREPSStatistics(csv, index);
                ep = ep + 1;
            end
            
        end
        
        clearvars csv
        
        
        J_history = [];
        for o = 1:length(data)
            J(o,count) = mean([data(o).policies.J]);
            %             J_history = [J_history, data(k).J'];
        end
        testname{count} = [algorithms{i}];
        count = count + 1;
    end
end
% hold off;
%%
figure(2);
hold on;
for i = 1:length(testname)
    plot(smooth(J(:,i)), 'Linewidth', 1.5)
    %     plot(J(:,i), 'Linewidth', 1.5)
    %     disp(testname{i})
    %     axis([0,nbUpdates,-400,-100])
    %     pause
end
grid on;
title([domain,', episodes:' num2str(nbEpisodes),', iter:',num2str(nbUpdates)]);
legend(testname, 'location', 'southeast');
set(gca,'yscale','log');
hold off;

%% Plot results
% figure; shadedErrorBar(1:size(J_history,2), ...
%     mean(J_history), ...
%     2*sqrt(diag(cov(J_history))), ...
%     {'LineWidth', 2'}, 1);
xlabel('Iterations')
ylabel('Average return')