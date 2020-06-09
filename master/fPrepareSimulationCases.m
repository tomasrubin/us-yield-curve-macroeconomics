function [simulation_cases,simulation_runs_iCase,simulation_runs_process,simulation_runs_response_process,simulation_runs_response_process_shape] = fPrepareSimulationCases()

%% set some constant concerning the simulations
nRandi_min = 0; % the minimal number of points per curve - the same for all cases
simulation_variables.nGridTime = [300 600 900 1200]; % maximum time horizon T
simulation_variables.nRandi_max = [10 20 40 60]; % maximum sampled points per curve


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  prepare the structure for different simulation cases

nCases = length(simulation_variables.nGridTime) * length(simulation_variables.nRandi_max);

simulation_cases = [];
simulation_cases.nGridTime = zeros(1,nCases);
simulation_cases.nRandi_max = zeros(1,nCases);
simulation_cases.memory = zeros(1,nCases);

ii = 1;
for i1 = simulation_variables.nRandi_max
    for i2 = simulation_variables.nGridTime    
            simulation_cases.nRandi_max(ii) = i1;
            simulation_cases.nGridTime(ii) = i2;            
            simulation_cases.memory(ii) = ((i1/10)^2 * i2/300) / 6 + 4;
            ii = ii + 1;
    end
end

% reshape( simulation_cases.nGridTime, 4, 4)
% reshape( simulation_cases.memory, 4, 4)




simulation_runs_iCase = {};
simulation_runs_process = {};
simulation_runs_response_process = {};
simulation_runs_response_process_shape = {};

for iCase = 1:nCases
    
    if simulation_cases.memory(iCase) >= 20
        
        % individually        
        for process = [0 1 4]
            for response = [1 2 3]
                for shape = [2 3]
                    simulation_runs_iCase{end+1} = iCase;
                    simulation_runs_process{end+1} = process;
                    simulation_runs_response_process{end+1} = response;
                    simulation_runs_response_process_shape{end+1} = shape; 
                end
            end
        end               
        
    elseif simulation_cases.memory(iCase) >= 10
        % both shapes together
        for process = [0 1 4]
            for response = [1 2 3]
                simulation_runs_iCase{end+1} = iCase;
                simulation_runs_process{end+1} = process;
                simulation_runs_response_process{end+1} = response;
                simulation_runs_response_process_shape{end+1} = [2 3];                
            end
        end
        
    elseif simulation_cases.memory(iCase) >= 6
        % all shapes and responses together
        for process = [0 1 4]
            simulation_runs_iCase{end+1} = iCase;
            simulation_runs_process{end+1} = process;
            simulation_runs_response_process{end+1} = [1 2 3];
            simulation_runs_response_process_shape{end+1} = [2 3];
        end
        
    else
        % all together
        simulation_runs_iCase{end+1} = iCase;
        simulation_runs_process{end+1} = [0 1 4];
        simulation_runs_response_process{end+1} = [1 2 3];
        simulation_runs_response_process_shape{end+1} = [2 3]; 
    end
    
    
end

nRuns = length(simulation_runs_iCase);