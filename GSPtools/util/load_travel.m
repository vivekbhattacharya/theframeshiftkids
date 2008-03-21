function [Travel] = load_travel()
    load Travel2.mat;
    Travel.uua = Travel.cuu * 4.7;
end