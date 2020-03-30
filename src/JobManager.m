classdef JobManager < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        jobCounter = 1;
        jobs = Job.empty;
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
            % Run all jobs
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
            
            % Create table for the options of the plots
%             fig = uifigure;
%             fig.Position(3:4) = [1400 360];
%             uit = uitable(fig);
%             uit.Position = [20 20 620 310];
%             tObject = struct2table([obj.jobs.problemConfiguration]);
%             tObject.material_1 = struct2table(tObject.material_1);
%             tObject.material_2 = struct2table(tObject.material_2);
%             uit.Data = tObject;
%             uit.RowName = {obj.jobs.name};
%             uit.ColumnWidth = {50, 40, 40, 60, 'auto', 'auto', 80};
%             
%             uit = uitable(fig);
%             uit.Position = [650 20 500 310];
%             asArray = length(obj.jobs) == 1;
%             tObject = struct2table([obj.jobs.femConfiguration], 'AsArray', asArray);
%             uit.Data = tObject;
            
            for job = obj.jobs
                job.plotResult();
            end
        end
        
        function folderPath = getFolderPath(obj)
            subFolderName = ['jobs', datestr(now,'yyyy-mm-dd HH_MM_SS-FFF'), '/'];
            mkdir(obj.savePath, subFolderName);
            folderPath = [obj.savePath, subFolderName];
        end
    end
end
