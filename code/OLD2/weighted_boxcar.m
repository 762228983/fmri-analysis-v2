function [weighted_boxcar_matrix, regressor_names,t] = weighted_boxcar(para_file, weight_file, sr)

[boxcar_matrix, condition_names, t] = boxcar_from_para(para_file, sr);


weighted_boxcar_matrix = boxcar_matrix * boxcar_weights;
regressor_names = W.regressor_names;