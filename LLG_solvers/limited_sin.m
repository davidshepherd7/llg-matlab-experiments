function y = limited_sin(x)
% As sin(x) but with a minimum value as specified by min_val.

% Example use: to prevent (meaningless) divergence in LLG as theta -> 0.

min_val = 0.001;
y = sin(x);

if y < min_val
    y = min_val;
end