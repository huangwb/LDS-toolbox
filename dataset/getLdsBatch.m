function lds = getLdsBatch(data,LDS_options)
%GETLDSBATCH Summary of this function goes here
%   Detailed explanation goes here

    function lds = cellDytex(data_i)
        lds = dytex(data_i, LDS_options);
    end

lds = cellfun(@cellDytex,data,'UniformOutput',false);

% lds = cell(1,blocks.number);
% for k=1:blocks.number
%     lds{k} = dytex(blocks.data{k}, LDS_options);
% end

end

