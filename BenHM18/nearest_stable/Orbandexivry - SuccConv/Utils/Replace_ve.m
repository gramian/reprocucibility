infilename  = 'testfile.txt'
outfilename = 'outtest.txt'

fid = fopen(infilename,'r')

if fid == -1
    error('Could not open the input file : %s',infilename)
end
fid2 = fopen(outfilename,'w')
if fid2 == -1
    fclose(fid);
    error('Could not open the output file : %s',outfilename)
end

try
    %read lines
    isline = 1;
    tline = fgetl(fid);
    if ~ischar(tline),  error('Problem reading file in replace_ve'),   end
    while isline
        %find begining of command
        K = strfind(tline,'\ve{');
        if ~isempty(K)
            indmid = zeros(size(K));
            indend = zeros(size(K));
            for k=1:length(K)
                %finding middle }{
                m = K(k)+4;
                countbrace = 0;
                notmiddle = 1;
                while notmiddle
                    if strcmp(tline(m),'{')
                        countbrace = countbrace+1;
                    elseif strcmp(tline(m),'}')
                        if countbrace == 0;
                            if strcmp(tline(m+1),'{')
                                notmiddle = 0;
                                indmid(k) = m;
                            end
                        else
                            countbrace = countbrace-1;
                        end
                        if countbrace < 0
                            error('Badcounting in the middle')
                        end
                    end
                    m = m +1;
                end
                %finding last }
                countbracee = 0;
                notend = 1;
                e = indmid(k)+2;
                while notend
                    if strcmp(tline(e),'{')
                        countbracee = countbracee + 1;
                    elseif strcmp(tline(e),'}')
                        if countbracee == 0
                            notend = 0;
                            indend(k) = e;
                        else
                            countbracee = countbracee -1;
                        end
                        if countbracee < 0
                            error('Badcounting at the end')
                        end
                    end
                    e = e + 1;
                end
            end
        end
        
        newline = tline;
        if ~isempty(K)
            
            %regroup indices
            indmat = [K;indmid;indend];
            
            %Construct newline
            for i = length(K):-1:1
                %Start from the end of the line
                
                %replace final }
                replend = '\rangle ';
                newline = [newline(1:(indmat(3,i)-1)) replend newline((indmat(3,i)+1):end)];
                                
                %replace middle }{
                replmid = ',';
                newline = [newline(1:(indmat(2,i)-1)) replmid newline((indmat(2,i)+2):end)];
                
                %replace beginning of \ve{
                replbeg = '\text{\red{TODO}} \langle ';
                newline = [newline(1:(indmat(1,i)-1)) replbeg newline((indmat(1,i)+4):end)];
            end
        end
        %prinf in file
        fprintf(fid2,'%s\n',newline);
        
        %Get a line in the original file
        tline = fgetl(fid);
        if ~ischar(tline), isline = 0; end
    end
    fclose(fid);
    fclose(fid2);
catch
    fclose(fid);
    fclose(fid2);
    rethrow(lasterror);
end
