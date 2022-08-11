function eyelink_close(edfFile)

    % a function that used when quit the presentation loop 
 
    % Stop recording and close file
    Eyelink('StopRecording');
    Eyelink('CloseFile');
    
    % Download edf data
    try 
        fprintf('Receiving data file ''%s''\n', edfFile);
        statusFlag = Eyelink('ReceiveFile');
        if statusFlag > 0
            fprintf('ReceiveFile status %d\n', statusFlag);
        end
        if 2 == exist(edfFile, 'file')
            fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, pwd);
        end
    catch 
        fprintf('Problem receiving data file ''%s''\n', edfFile );
    end
    
    eyelinkdata = strcat(edfFile,'.edf');
    eyepath     = 'eyelinkdata';
    copyfile(eyelinkdata,eyepath);
    
    
    % Shut down the eyetracker
    Eyelink('Shutdown');


end