function [colors]=define_colors(q,LL)
if any(find(q==LL))
    colors = [0 0 0; % black for -1
        1 0 0; % red for 0
        1 0.5 0; % orange for 1
        1 1 0; % yellow for 2
        0 1 0; % green for 3, etc.
        0.5 0 1; %purple for 4
        0 0.2 1; % blue for 5
        0.2 0.8 1; % cyan for 6
        1 0 1];% pink for 7
else
    colors = [%0 0 0; % black for -1
        1 0 0; % red for 0
        1 0.5 0; % orange for 1
        1 1 0; % yellow for 2
        0 1 0; % green for 3, etc.
        0.5 0 1; %purple for 4
        0 0.2 1; % blue for 5
        0.2 0.8 1; % cyan for 6
        1 0 1];% pink for 7
end