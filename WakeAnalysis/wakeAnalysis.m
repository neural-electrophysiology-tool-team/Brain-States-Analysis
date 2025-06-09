

classdef wakeAnalysis < sleepAnalysis

    methods

        function locTable = getHeadPosition(obj)
            % read the parquet 
            videoPath = obj.recTable.VideoFiles(obj.currentPRec);
            [videosFolderPath,vidName,~] = fileparts(videoPath{1});
            locFilename = strcat(videosFolderPath,'/predictions/front_head_ephys_resnet_101__',vidName,'.parquet');

            % check file exists:
            if ~isfile(locFilename)
                error('File does not exist: %s', locFilename);
            end
            locTable = parquetread(locFilename);
            ArenaCSV = obj.getArenaCSVs();
            locTable.t_ms = ArenaCSV.oeCamTrigs(1:height(locTable));


            obj.files.ArenaLocation = [obj.currentAnalysisFolder filesep 'arenaLocation.mat'];
            save(obj.files.ArenaLocation,"locTable",'-mat')

        end

    end

end