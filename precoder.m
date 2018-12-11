classdef precoder
    methods 
        function obj = precoder(varargin)
            p = inputParser;

            p.KeepUnmatched = true;

            p.parse(varargin{:});
            params = p.Results;
        end

        function [rate] = simulate(obj, h, S, P, N0, MAI, iterations)
            rate = zeros(iterations, 1);
        end
    end
end
