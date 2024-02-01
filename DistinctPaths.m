function [PATH] = DistinctPaths(n,s,sp)
% compute matrix containing all the distinct paths for a model with s seasons
% in the annual cycle and n habitats, where seasons indexed by the entries
% of sp have paths specified. 
% PATH lists the habitats that the focal subpopulation travels to. If an
% entry in PATH equals 0, then the habitat is unspecified and so the
% population travelling to all habitats during the relevant season will be
% tracked. 

%% PATH IN ONE PARTICULAR SEASON IS SPECIFIED 
if length(sp) == 1
    % POPULATE MATRIX WITH ALL UNIQUE PATHS
    count = 1;
    PATH = zeros(n^(length(sp)+1),s+1); 
    for jj = 1:n
        for ii = 1:n 
            PATH(count,sp) = jj; 
            PATH(count,sp+1) = ii; 
            count = count +1;
        end
    end
%% TWO SEASONS CONTAIN A SPECIFIED PATH     
elseif length(sp) == 2
    % POPULATE MATRIX WITH ALL UNIQUE PATHS
    count = 1; 
    % If the two seasons with specified paths ARE consecutive
    if ismember(sp(1)+1,sp) == 1 
        % initiate matrix to store all unique paths
        PATH = zeros(n^(length(sp)+1),s+1); 
        % populate PATH with all possible paths
        for jj = 1:n 
            for ii = 1:n 
                for hh = 1:n 
                    PATH(count,sp(1)) = jj; 
                    PATH(count,sp(1)+1) = ii;
                    PATH(count,sp(1)+2) = hh;
                    count = count+1; 
                end
            end
        end
    % If the two seasons with specified paths are NOT consecutive    
    else %if ismember(sp(1)+1,sp) == 0
        % initiate matrix to store all unique paths
        PATH = zeros(n^(length(sp)+2),s+1); 
        % populate PATH with all possible paths
        for jj = 1:n 
            for ii = 1:n
                for hh = 1:n
                    for gg = 1:n 
                        PATH(count,sp(1)) = jj; 
                        PATH(count,sp(1)+1) = ii; 
                        %
                        PATH(count,sp(2)) = hh; 
                        PATH(count,sp(2)+1) = gg; 
                        count = count+1;
                    end
                end
            end
        end        
    end
%% THREE SEASONS CONTAIN A SPECIFIED PATH
elseif length(sp) == 3
    % POPULATE A MATRIX WITH ALL UNIQUE PATHS
    count = 1; 
    % If the three seasons with specified paths ARE all consecutive
    if ismember(sp(1)+1,sp) == 1 && ismember(sp(1)+2,sp) == 1
        % initiate matrix to store all unique paths
        PATH = zeros(n^(length(sp)+1),s+1); 
        % populate PATH with all possible paths
        for jj = 1:n 
            for ii = 1:n 
                for hh = 1:n 
                    for gg = 1:n
                        PATH(count,sp(1)) = jj; 
                        PATH(count,sp(1)+1) = ii;
                        PATH(count,sp(1)+2) = hh;
                        PATH(count,sp(1)+3) = gg; 
                        count = count+1; 
                    end
                end
            end
        end
    % If ONLY the first two seasons with specified paths ARE consecutive,
    % requires entries of sp to be specified in increasing order
    elseif ismember(sp(1)+1,sp) == 1 
        % initiate matrix to store all unique paths
        PATH = zeros(n^(length(sp)+2),s+1); 
        % populate PATH with all possible paths
        for jj = 1:n 
            for ii = 1:n
                for hh = 1:n
                    for gg = 1:n 
                        for ff = 1:n
                            PATH(count,sp(1)) = jj; 
                            PATH(count,sp(1)+1) = ii; 
                            PATH(count,sp(1)+2) = hh; 
                            %
                            PATH(count,sp(3)) = gg; 
                            PATH(count,sp(3)+1) = ff; 
                            count = count+1;
                        end
                    end
                end
            end
        end
    % If ONLY the last two seasons with specified paths ARE consecutive,
    % requires entries of sp to be specified in increasing order
    elseif ismember(sp(2)+1,sp) == 1
        % initiate matrix to store all unique paths
        PATH = zeros(n^(length(sp)+2),s+1); 
        % populate PATH with all possible paths
        for jj = 1:n 
            for ii = 1:n
                for hh = 1:n
                    for gg = 1:n 
                        for ff = 1:n
                            PATH(count,sp(1)) = jj; 
                            PATH(count,sp(1)+1) = ii; 
                            %
                            PATH(count,sp(2)) = hh; 
                            PATH(count,sp(2)+1) = gg; 
                            PATH(count,sp(2)+2) = ff; 
                            count = count+1;
                        end
                    end
                end
            end
        end
    % If all the seasons with specified paths are NOT consecutive
    else
        % initiate matrix to store all unique paths
        PATH = zeros(n^(length(sp)+3),s+1); 
        % populate PATH with all possible paths
        for jj = 1:n 
            for ii = 1:n
                for hh = 1:n
                    for gg = 1:n 
                        for ff = 1:n
                            for ee = 1:n 
                                PATH(count,sp(1)) = jj; 
                                PATH(count,sp(1)+1) = ii; 
                                %
                                PATH(count,sp(2)) = hh; 
                                PATH(count,sp(2)+1) = gg;
                                %
                                PATH(count,sp(3)) = ff;
                                PATH(count,sp(3)+1) = ee; 
                                count = count+1;
                            end
                        end
                    end
                end
            end
        end
    end
%% FOUR SEASONS CONTAIN A SPECIFIED PATH
elseif length(sp) == 4
   % POPULATE MATRIX WITH ALL UNIQUE PATHS
   count = 1;
   % If the four seasons with specified paths ARE all consecutive
   if ismember(sp(1)+1,sp) == 1 && ismember(sp(1)+2,sp) == 1 && ismember(sp(1)+3,sp) == 1
       % initiate matrix to store all unique paths 
       PATH = zeros(n^(length(sp)+1),s+1);
       % populate PATH with all possible paths
       for jj = 1:n 
            for ii = 1:n 
                for hh = 1:n 
                    for gg = 1:n
                        for ff = 1:n 
                            PATH(count,sp(1)) = jj; 
                            PATH(count,sp(1)+1) = ii;
                            PATH(count,sp(1)+2) = hh;
                            PATH(count,sp(1)+3) = gg; 
                            PATH(count,sp(1)+4) = ff; 
                            count = count+1; 
                        end 
                    end
                end
            end
        end
   % If ONLY the first three seasons with specified paths ARE consecutive,
   % requires entries of sp to be specified in increasing order
   elseif ismember(sp(1)+1,sp) == 1 && ismember(sp(1)+2,sp) == 1
       % initaite matrix to store all unique paths
       PATH = zeros(n^(length(sp)+2),s+1); 
       % populate PATH with all possible paths
       for jj = 1:n 
            for ii = 1:n
                for hh = 1:n
                    for gg = 1:n 
                        for ff = 1:n
                            for ee = 1:n
                                PATH(count,sp(1)) = jj; 
                                PATH(count,sp(1)+1) = ii;
                                PATH(count,sp(1)+2) = hh;
                                PATH(count,sp(1)+3) = gg; 
                                %
                                PATH(count,sp(4)) = ff;
                                PATH(count,sp(4)+1) = ee; 
                                count = count+1;
                            end
                        end
                    end
                end
            end
       end
   % If ONLY the last three seasons with specified paths ARE consecutive,
   % requires entries of sp to be specified in increasing order
   elseif ismember(sp(2)+1,sp) == 1 && ismember(sp(2)+2,sp) == 1
       % initaite matrix to store all unique paths
       PATH = zeros(n^(length(sp)+2),s+1); 
       % populate PATH with all possible paths
       for jj = 1:n 
            for ii = 1:n
                for hh = 1:n
                    for gg = 1:n 
                        for ff = 1:n
                            for ee = 1:n
                                PATH(count,sp(1)) = jj; 
                                PATH(count,sp(1)+1) = ii;
                                %
                                PATH(count,sp(2)) = hh;
                                PATH(count,sp(2)+1) = gg; 
                                PATH(count,sp(2)+2) = ff;
                                PATH(count,sp(2)+3) = ee; 
                                count = count+1;
                            end
                        end
                    end
                end
            end
       end
   % If ONLY the first two seasons with specified paths ARE consecutive,
   % requires entries of sp to be specified in increasing order
   elseif ismember(sp(1)+1,sp) == 1 && ismember(sp(2)+1,sp) == 0 && ismember(sp(3)+1,sp) == 0
       % initaite matrix to store all unique paths
       PATH = zeros(n^(length(sp)+3),s+1); 
       % populate PATH with all possible paths
       for jj = 1:n 
            for ii = 1:n
                for hh = 1:n
                    for gg = 1:n 
                        for ff = 1:n
                            for ee = 1:n
                                for dd = 1:n 
                                    PATH(count,sp(1)) = jj; 
                                    PATH(count,sp(1)+1) = ii;
                                    PATH(count,sp(1)+2) = hh;
                                    %
                                    PATH(count,sp(3)) = gg; 
                                    PATH(count,sp(3)+1) = ff;
                                    %
                                    PATH(count,sp(4)) = ee; 
                                    PATH(count,sp(4)+1) = dd; 
                                    count = count+1;
                                end
                            end
                        end
                    end
                end
            end
       end
    % If the first two seasons AND the last two seasons with specified 
    % paths ARE consecutive, requires entries of sp to be specified in 
    % increasing order
    elseif ismember(sp(1)+1,sp) == 1 && ismember(sp(2)+1,sp) == 0 && ismember(sp(3)+1,sp) == 1
        % initiate matrix to store all unique paths
        PATH = zeros(n^(length(sp)+2),s+1); 
        % populate PATH with all possible paths
        for jj = 1:n 
            for ii = 1:n
                for hh = 1:n
                    for gg = 1:n 
                        for ff = 1:n
                            PATH(count,sp(1)) = jj; 
                            PATH(count,sp(1)+1) = ii;
                            PATH(count,sp(1)+2) = hh; 
                            %
                            PATH(count,sp(2)) = hh;
                            PATH(count,sp(2)+1) = gg; 
                            PATH(count,sp(2)+2) = ff;  
                            count = count+1;
                        end
                    end
                end
            end
        end
   % If ONLY the mid two seasons with specified paths ARE consecutive,
   % requires entries of sp to be specified in increasing order
   elseif ismember(sp(1)+1,sp) == 1 && ismember(sp(2)+1,sp) == 0 && ismember(sp(3)+1,sp) == 0
       % initaite matrix to store all unique paths
       PATH = zeros(n^(length(sp)+3),s+1); 
       % populate PATH with all possible paths
       for jj = 1:n 
            for ii = 1:n
                for hh = 1:n
                    for gg = 1:n 
                        for ff = 1:n
                            for ee = 1:n
                                for dd = 1:n 
                                    PATH(count,sp(1)) = jj; 
                                    PATH(count,sp(1)+1) = ii;
                                    %
                                    PATH(count,sp(2)) = hh;
                                    PATH(count,sp(2)+1) = gg; 
                                    PATH(count,sp(2)+2) = ff;
                                    %
                                    PATH(count,sp(4)) = ee; 
                                    PATH(count,sp(4)+1) = dd; 
                                    count = count+1;
                                end
                            end
                        end
                    end
                end
            end
       end 
   % If ONLY the last two seasons with specified paths ARE consecutive,
   % requires entries of sp to be specified in increasing order
   elseif ismember(sp(1)+1,sp) == 0 && ismember(sp(2)+1,sp) == 0 && ismember(sp(3)+1,sp) == 1
       % initaite matrix to store all unique paths
       PATH = zeros(n^(length(sp)+3),s+1); 
       % populate PATH with all possible paths
       for jj = 1:n 
            for ii = 1:n
                for hh = 1:n
                    for gg = 1:n 
                        for ff = 1:n
                            for ee = 1:n
                                for dd = 1:n 
                                    PATH(count,sp(1)) = jj; 
                                    PATH(count,sp(1)+1) = ii;
                                    %
                                    PATH(count,sp(2)) = hh; 
                                    PATH(count,sp(2)+1) = gg;
                                    %
                                    PATH(count,sp(3)) = ff; 
                                    PATH(count,sp(3)+1) = ee; 
                                    PATH(count,sp(3)+2) = dd;
                                    count = count+1;
                                end
                            end
                        end
                    end
                end
            end
       end
    % If all seasons with a specified path are NOT consecutive
   elseif ismember(sp(1)+1,sp) == 0 && ismember(sp(2)+1,sp) == 0 && ismember(sp(3)+1,sp) == 0 
       PATH = zeros(n^(length(sp)+4),s+1); 
       % populate PATH with all possible paths
       for jj = 1:n 
            for ii = 1:n
                for hh = 1:n
                    for gg = 1:n 
                        for ff = 1:n
                            for ee = 1:n
                                for dd = 1:n 
                                    for cc = 1:n 
                                        PATH(count,sp(1)) = jj; 
                                        PATH(count,sp(1)+1) = ii;
                                        %
                                        PATH(count,sp(2)) = hh; 
                                        PATH(count,sp(2)+1) = gg;
                                        %
                                        PATH(count,sp(3)) = ff; 
                                        PATH(count,sp(3)+1) = ee; 
                                        % 
                                        PATH(count,sp(4)) = dd;
                                        PATH(count, sp(4)+1) = cc; 
                                        count = count+1;
                                    end
                                end
                            end
                        end
                    end
                end
            end
       end
   end
%% TWELVE SEASONS CONTAIN A SPECIFIED PATH
elseif length(sp) == 12
   % POPULATE MATRIX WITH ALL UNIQUE PATHS
   count = 1;
   % For our purposes, the most seasons the annual cycle is split up into
   % is  12, so we need only consider the case where ... 
   % ALL seasons with specified paths ARE consecutive
       % initiate matrix to store all unique paths 
       PATH = zeros(n^(length(sp)+1),s+1);
       % populate PATH with all possible paths
       for jj = 1:n 
        for ii = 1:n 
         for hh = 1:n 
          for gg = 1:n
           for ff = 1:n 
            for ee = 1:n 
             for dd = 1:n 
              for cc = 1:n 
               for bb = 1:n
                for aa = 1:n
                 for zz = 1:n 
                  for yy = 1:n
                   for xx = 1:n 
                       PATH(count,sp(1)) = jj; 
                       PATH(count,sp(1)+1) = ii;
                       PATH(count,sp(1)+2) = hh;
                       PATH(count,sp(1)+3) = gg; 
                       PATH(count,sp(1)+4) = ff; 
                       PATH(count,sp(1)+5) = ee; 
                       PATH(count,sp(1)+6) = dd;
                       PATH(count,sp(1)+7) = cc;
                       PATH(count,sp(1)+8) = bb; 
                       PATH(count,sp(1)+9) = aa;
                       PATH(count,sp(1)+10) = zz; 
                       PATH(count,sp(1)+11) = yy;
                       PATH(count,sp(1)+12) = xx;
                       count = count+1; 
                   end
                  end
                 end
                end
               end
              end
             end
            end
           end
          end
         end
        end
       end
end
end