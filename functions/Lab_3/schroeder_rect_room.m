function fs = schroeder_rect_room(room,TR60,varargin)
    % Compute the schroeder frequency for a given rectangular room
    %
    % Parameters:
    % room - Dimensions of the rectangular room [lx,ly,lz]
    % RT60 - Averaged reverberation time acrross frequencies
    % Optional parameters:
    %   'c0' - sound speed in the medium
    %
    % Returns:
    % fs - schroeder frequency, (modal density >3)

    % -------------------------------- init function -----------------------------
    p = inputParser;
    addRequired(p, 'room');
    addRequired(p, 'TR60');
    addParameter(p, 'c0', 343, @isnumeric)

    parse(p, room, TR60, varargin{:})

    c0 = p.Results.c0;

    V = prod(room);

    fs = sqrt(c0^3 * TR60 / (4 * log(10) * V));
end