function output_prb=repeat_prb(Ka,M)

no_collision=exp(sum(log(1-(0:Ka-1)./M)));
one_doublecollision=nchoosek(Ka,2)*(1/M)*exp(sum(log(1-(1:Ka-2)./M)));
triple_collision=nchoosek(Ka,3)*(1/M^2)*exp(sum(log(1-(1:Ka-3)./M)));
two_doublecollision=0.5*nchoosek(Ka,2)*nchoosek(Ka-2,2)*(1/M^2)*(1-1/M)*exp(sum(log(1-(2:Ka-3)./M))); %this is approximate


output_prb=[no_collision one_doublecollision two_doublecollision triple_collision];

end

