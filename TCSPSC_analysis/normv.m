function [normv] = normalise(v)

normv = (v-min(v))./(max(v)-min(v));