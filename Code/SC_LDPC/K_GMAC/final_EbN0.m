function output=final_EbN0(Eb_N01,Eb_N02,B1,B2)

output=10*log10((B1*10^(Eb_N01/10)+B2*10^(Eb_N02/10))/(B1+B2));

end
