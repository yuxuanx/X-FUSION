function [ x_est ] = state_extract2(r_update, x_update, existThresh )
%Best state extraction

idices = (r_update > existThresh);
x_est = x_update(:,idices);

end

