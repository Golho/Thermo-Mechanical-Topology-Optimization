classdef JobManager < handle
    %JOBMANAGER Manages the TO jobs
    
    properties
        jobCounter = 1;         
        jobs = Job.empty;
        %logger;
        savePath = "results/";
    end
    
    methods
        function obj = JobManager()
            %obj.logger = logger;
        end
        
        function add(obj, job)
            obj.jobs(end+1) = job;
            
            % If the job doesn't have a name, assign it default name with
            % counting suffix
            if isempty(job.name)
                job.name = ['job', num2str(obj.jobCounter)];
                obj.jobCounter = obj.jobCounter + 1;
            end
        end
        
        function loadAll(obj, folderName)
            if nargin <= 1
                folderName = pwd;
            end
            
            matFiles = dir(fullfile(folderName, '*.mat'));
            for matFile = matFiles'
                s = load(fullfile(matFile.folder, matFile.name));
                if isa(s.obj, 'Job')
                    obj.add(s.obj);
                end
            end
        end
        
        function runAndSaveAll(obj)
            folderPath = obj.getFolderPath();
            for job = obj.jobs
                job.run();
                job.save(folderPath);
            end
        end
        
        function runAll(obj)
            for job = obj.jobs
                job.run();
            end
        end
        
        function saveAll(obj)
            folderPath = obj.getFolderPath();
            for job = obj.jobs
                job.save(folderPath)
            end
        end
        
        function plotAll(obj)
            for job = obj.jobs
                job.plotResult();
            end
        end
    end
    
    methods(Access = protected)
        function folderPath = getFolderPath(obj)
            subFolderName = "jobs" + datestr(now,'yyyy-mm-dd HH_MM_SS-FFF') + "/";
            mkdir(obj.savePath, subFolderName);
            folderPath = obj.savePath + subFolderName;
        end
    end
end
