function [codes] = get_color(varargin)
% get_color input is a request for color for some category
% output is the RGB color code
    if length(varargin) > 1
        codes = zeros(length(varargin),3);
        names = [varargin{:}];
    else
        codes = zeros(length(varargin{1}),3);
        names = varargin{1};
    end
    % Codes
    cds = {"Enhanced", [1 0.1529 0.1725]
            "Suppressed", [0.0157 0.0980 1]
            "None",[0.2627 0.8 0.2353]
            "All", [0.3 0.3 0.3]
            "Hit", [0.4784 0.6510 0.6157]
            "Miss", [0.7608 0.5843 0.7216]
            "FA", [0.7098 0.6588 0.2353]
            "mix Enhanced",[0.8706 0.4941 0.5373]
            "mix Suppressed", [0.6627 0.6431 0.9804]
            "mix None", [0.4196 0.7686 0.3843]};
        
    % return code
    for i = 1:length(names)
        row_ind = find(strcmpi([cds{:,1}],names(i)));
        codes(i,:) = cds{row_ind,2};
    end
end