function [qam8, qam_binary_map] = generate_qam8_map()
%GENERATE_QAM8_MAP Generate a rectangular 8-QAM constellation and bit map
%
% Outputs:
%   qam8           : 8x1 complex constellation points
%   qam_binary_map : 8x3 binary labels for each symbol
%
% Constellation shape:
%   I levels: [-3 -1 1 3]
%   Q levels: [-1 1]
%
% Gray-like labeling:
%   I-axis uses 2 bits: 00, 01, 11, 10
%   Q-axis uses 1 bit : 0, 1

    % Amplitude levels
    I_levels = [-3 -1 1 3];
    Q_levels = [-1 1];

    % Gray code for I-axis (2 bits)
    I_bits = [
        0 0
        0 1
        1 1
        1 0
    ];

    % Gray code for Q-axis (1 bit)
    Q_bits = [
        0
        1
    ];

    qam8 = zeros(8,1);
    qam_binary_map = zeros(8,3);

    idx = 1;
    for q = 1:length(Q_levels)
        for i = 1:length(I_levels)
            qam8(idx) = I_levels(i) + 1i*Q_levels(q);
            qam_binary_map(idx,:) = [I_bits(i,:) Q_bits(q)];
            idx = idx + 1;
        end
    end
end