function [y]=BinaryEntropy(p)
y=-(p.*log2(p)+(1-p).*log2(1-p));