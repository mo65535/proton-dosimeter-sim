function [ insiders, starters, stoppers, crossers, no_shows ] = load_mc_data( file_path )
%LOAD_MC_DATA Load Monte Carlo data from text file

    fp = fopen(file_path, 'r');

    C = textscan(fp, '%s', 'delimiter','\n');

    insiders = [];
    starters = [];
    stoppers = [];
    crossers = [];
    no_shows = [];

    N_insiders = -1;
    N_starters = -1;
    N_stoppers = -1;
    N_crossers = -1;
    N_no_shows = -1;

    i = 1;
    while (N_insiders == -1) || ...
          (N_starters == -1) || ...
          (N_stoppers == -1) || ...
          (N_crossers == -1) || ...
          (N_no_shows == -1)

        [match, match_count] = sscanf(C{1}{i}, ' Number of  %s:%d');

        if match_count == 1
          % Trim ':' off match
          match = match(1:end-1);

          % get count
          N = sscanf(C{1}{i}, sprintf(' Number of  %s: %%d',match));

          %fprintf('match = %s, N = %d\n', match, N);

          if strcmp('insiders',match)
              N_insiders = N;
          end
          if strcmp('starters',match)
              N_starters = N;
          end
          if strcmp('stoppers',match)
              N_stoppers = N;
          end
          if strcmp('crossers',match)
              N_crossers = N;
          end
          if strcmp('no-shows',match)
              N_no_shows = N;
          end
        end

        i = i+1;
    end

    data_found = false;

    while data_found == false
        if strcmp('E_dep_B by insiders [eV]:',C{1}{i})
            data_found = true;
        end
        i = i+1;
    end


    reformat = cellfun(@str2double, C{1}(i:i+N_insiders-1));
    if size(reformat) == [0 1]
        reformat = [];
    end
    insiders = reformat;
    i = i+N_insiders+1;


    reformat = cellfun(@str2double, C{1}(i:i+N_starters-1));
    if size(reformat) == [0 1]
        reformat = [];
    end
    starters = reformat;
    i = i+N_starters+1;


    reformat = cellfun(@str2double, C{1}(i:i+N_stoppers-1));
    if size(reformat) == [0 1]
        reformat = [];
    end
    stoppers = reformat;
    i = i+N_stoppers+1;


    reformat = cellfun(@str2double, C{1}(i:i+N_crossers-1));
    if size(reformat) == [0 1]
        reformat = [];
    end
    crossers = reformat;
    i = i+N_crossers+1;


    reformat = cellfun(@str2double, C{1}(i:i+N_no_shows-1));
    if size(reformat) == [0 1]
        reformat = [];
    end
    no_shows = reformat;
    i = i+N_no_shows+1;

end

