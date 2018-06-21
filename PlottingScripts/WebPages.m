% Change to the appropriate directory
allPathname = uigetdir('D:/Imaging','Select the directory');
cd(allPathname);
findDirs = dir;

% Make the main page
fileID = fopen(strcat(allPathname,'\Overview','.htm'),'w');
fprintf(fileID,'<html>\n');
fprintf(fileID,'<head>\n<meta http-equiv="Content-Type" content="text/html; charset=windows-1252">\n');
fprintf(fileID,'<title>%s</title>\n','Overview');
fprintf(fileID,'<head>\n');
fprintf(fileID,'<body topmargin="0" leftmargin="0" rightmargin="0" bottommargin="0" marginheight="0" marginwidth="0" bgcolor="#FFFFFF">\n');
for fID = 3:length(findDirs)
    fName = findDirs(fID).name; 
    if exist(fName) == 7
        subDir = dir(fName);
        for fID2 = 3:length(subDir)
            fName2 = subDir(fID2).name;
            if strcmp(fName2(end-7:end),'Traj.jpg')
                fprintf(fileID,'<p style="margin-left: 20"><b><font color="#B8C0F0" face="Arial" size="4">&nbsp;</font></b><font face="Arial" size="%f" color="#000000"> %s</font></p>',16,strcat(fName,', ',fName2));
                fprintf(fileID,'\n');
                strcat(fName,', ',fName2)
                expNotes = input('Notes for the experiment?: ');
                fprintf(fileID,'<p style="margin-left: 20"><b><font color="#B8C0F0" face="Arial" size="4">&nbsp;</font></b><font face="Arial" size="%f" color="#000000"> %s</font></p>',12,expNotes);
                fprintf(fileID,'\n');
                fprintf(fileID,'<A HREF="%s">\n',strcat(fName,'\',fName2(1:end-9),'.htm'));
                fprintf(fileID,'<IMG SRC="%s" ALT="%s" WIDTH=1200 HEIGHT=900>\n',strcat(fName,'\',fName2),strcat(fName,'\',fName2));
                fprintf(fileID,'</A>\n');
            end
        end
    end
end
fprintf(fileID,'</body>\n');
fprintf(fileID,'</html>\n');
fclose(fileID);

% Make the sub pages
for fID = 3:length(findDirs)
    fName = findDirs(fID).name; 
    if exist(fName) == 7
        subDir = dir(fName);
        for fID2 = 3:length(subDir)
            fName2 = subDir(fID2).name;
            if strcmp(fName2(end-7:end),'Traj.jpg')
                fileID = fopen(strcat(allPathname,'\',fName,'\',fName2(1:end-9),'.htm'),'w');
                fprintf(fileID,'<html>\n');
                fprintf(fileID,'<head>\n<meta http-equiv="Content-Type" content="text/html; charset=windows-1252">\n');
                fprintf(fileID,'<title>%s</title>\n','Overview');
                fprintf(fileID,'<head>\n');
                fprintf(fileID,'<body topmargin="0" leftmargin="0" rightmargin="0" bottommargin="0" marginheight="0" marginwidth="0" bgcolor="#FFFFFF">\n');
                for fID3 = 3:length(subDir)
                    fName3 = subDir(fID3).name;
                    if (strcmp(fName2(1:5),fName3(1:5)) & (strcmp(fName3(end-11:end),'Overview.jpg') || strcmp(fName3(end-11:end),'DiffComp.jpg') || strcmp(fName3(end-9:end),'CompBD.jpg')))
                        fprintf(fileID,'<IMG SRC="%s" ALT="%s" WIDTH=600 HEIGHT=450>\n',fName3,fName3);
                        fprintf(fileID,'</A>');
                    end
                end
                fprintf(fileID,'</body>\n');
                fprintf(fileID,'</html>\n');
                fclose(fileID);
            end
        end
    end
end