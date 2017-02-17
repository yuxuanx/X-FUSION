function [ w_update ] = normalize( w_update )

w_update = w_update/(sum(w_update)+eps);

end

