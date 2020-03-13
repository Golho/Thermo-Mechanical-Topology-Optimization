classdef JobManager < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        jobCounter = 1;
        jobs = {};
        logger;
        savePath = 'results/';
    end
    
    methods
        function obj = JobManager()
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            %obj.logger = logger;
        end
        
        function add(obj, job)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj.jobs{end+1} = job;
            
            % If the job doesn't have a name, assign it default name with
            % counting suffix
            if isempty(job.name)
                job.name = ['job', num2str(obj.jobCounter)];
                obj.jobCounter = obj.jobCounter + 1;
            end
        end
        
        function runAll(obj)
            % Run all jobs
            for i = 1:length(obj.jobs)
                obj.jobs{i}.run();
            end
        end
        
        function saveAll(obj)
            subFolderName = ['jobs', datestr(now,'yyyy-mm-dd HH_MM_SS-FFF'), '/'];
            mkdir(obj.savePath, subFolderName);
            folderPath = [obj.savePath, subFolderName];
            for i = 1:length(obj.jobs)
                obj.jobs{i}.save(folderPath)
            end
        end
        
        function plotAll(obj)
            for i = 1:length(obj.jobs)
                obj.jobs{i}.plotResult();
            end
        end
    end
end
